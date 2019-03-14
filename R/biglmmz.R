#------------------------
# Main function biglmmz
#------------------------

#' Low-rank LMM with a single random effect and residual random effect.
#'
#' Model definition:
#' y_i =  X_i beta + u_i + e_i
#'
#' u ~ N(0, s2*h2*G)
#' e ~ N(0, s2 I)
#'
#' var(y) = V = s2 * (h2*G + I)
#'
#' @import bigstatsr
#'
#' @export
biglmmz <- function(y, X, Z, scale = FALSE,
  REML = TRUE,
  compute_K = TRUE,
  store_mat = FALSE, 
  verbose = 0)
{
  mc <- match.call()

  if(verbose) { cat("biglmmz:\n") }

  ### args
  stopifnot(!missing(y))
  if(!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  if(missing(X)) {
    X <- matrix(1, nrow = length(y), ncol = 1)
  }

  stopifnot(!missing(Z))

  stopifnot(nrow(y) == nrow(Z))

  ### process Z -> get a scaled 'double' FBM
  Z0 <- Z
  Z <- FBM(nrow(Z0), ncol(Z0))
  if (verbose) cat(" - convert and scale Z if necessary\n")
  big_apply(Z, function(Z, ind) {
    if (scale) {
      # scale Z such a way that ZZ' = GRM
      Z0_part <- Z0[, ind]
      
      col_means <- colMeans(Z0_part)
      col_freq <- col_means / 2  # col_means = 2 * col_freq
      col_sd <- sqrt(2 * ncol(Z) * col_freq * (1 - col_freq))
      
      Z[, ind] <- sweep(sweep(Z0_part, 2, col_means, "-"), 2, col_sd , "/")
    } else {
      Z[, ind] <- Z0[, ind]
    }
    NULL
  })

  ### pre-compute K = Z'Z
  K <- NULL
  if(compute_K) {
    K <- crossprod(Z)
  }
  
  ### optimize
  if (verbose) cat(" - optimize\n")
  out <- optimize(biglr_ll_inline, c(0, 1),
    y = y, X = X, Z = Z, K = K, REML = REML, verbose = verbose,
    maximum = TRUE)

  r2 <- out$maximum
  ll <- out$objective
  convergence <- NA

  ### estimates of fixed effects
  if (verbose) cat(" - fixed effects\n")
  
  est <- biglr_fixef(r2, y, X, Z, K = K, REML = REML)
  
  coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  coef <- within(coef, z <- beta / se)
  
  ### estimates of random effects
  if (verbose) cat(" - random effects\n")
  
  gamma <- r2
  s2 <- est$s2
  comp <- s2 * c(gamma, 1 - gamma)
  
  ### return
  mod <- list(
    gamma = gamma, s2 = s2, comp = comp,
    est = est, coef = coef,
    REML = REML, store_mat = store_mat)

  if(store_mat) {
    mod <- c(mod, list(y = y, X = X, Z = Z))
  }

  mod$lmm <- list(r2 = r2, ll = ll, convergence = convergence, REML = REML)

  return(mod)
}

#-------------------------------
# Estimates of fixed effects
#-------------------------------

# The codes to be revised & optimized.
biglr_fixef <- function(gamma, y, Xmat, Z, s2, K = NULL, REML = TRUE)
{
  missing_s2 <- missing(s2)

  n <- length(y)
  k <- ncol(Xmat)
  
  nk <- `if`(REML, n - k, n)

  # pre-compute `K = Z'Z`
  if(is.null(K)) {
    K <- crossprod(Z)
  }

  # need to compute `s2`
  # then we can get beta/se taking into account `s2`: V = s2 (gamma ZZ' + (1-gamma) I)
  if(missing_s2) {
    comp <- c(gamma, 1 - gamma)
  
    # compute `XVX` & coeff. `b`, where `XVX` = X' Vi X
    XVt <- biglr_cprodMatInv(comp, Z, Xmat, K, transpose = TRUE) # Vi' X
    XVX <- crossprod(XVt, Xmat) 
    b <- solve(XVX, crossprod(XVt, y)) # (XVX)_inv XV y
  
    Rmat <- y - Xmat %*% b # r <- as.numeric(Rmat)
    rVt <- biglr_cprodMatInv(comp, Z, Rmat, K, transpose = TRUE) # Vi' r
    yPy <- crossprod(rVt, Rmat) 
    s2 <- as.numeric(yPy / nk)
  }
  
  # update `comp` from unscaled to scaled by `s2`
  comp <- s2 * c(gamma, 1 - gamma)
  
  XVt <- biglr_cprodMatInv(comp, Z, Xmat, K, transpose = TRUE) # Vi' X
  XVX <- crossprod(XVt, Xmat) 
  b <- solve(XVX, crossprod(XVt, y))
  bcov <- solve(XVX)
  
  ### return
  out <- list(s2 = s2, b = as.numeric(b), bcov = bcov)
}

#--------------------------------------------
# Log-likelihood function for low-rank LMM
#--------------------------------------------

biglr_ll_inline <- function(gamma, y, Xmat, Z, K = NULL, REML = TRUE, verbose = 0)
{
  if(!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  
  n <- length(y)
  k <- ncol(Xmat)

  nk <- `if`(REML, n - k, n)

  comp <- c(gamma, 1 - gamma)

  # pre-compute K = Z'Z
  if(is.null(K)) {
    K <- crossprod(Z)
  }
  
  # compute `Sigma_det_log` = det(V), where V = comp[1] ZZ' + comp[2] I
  diag(K) <- diag(K) + comp[2] / comp[1]
  log_det_K <- determinant(K, log = TRUE)
  Sigma_det_log <- as.numeric(log_det_K$modulus) + n*log(comp[2]) +
    ncol(Z)*log(comp[1] / comp[2])

  # compute `XVX` & coeff. `b`, where `XVX` = X' Vi X
  ZtX <- crossprod(Z, Xmat)
  XVX <- (crossprod(Xmat) - crossprod(ZtX, solve(K, ZtX)))

  XVty <- crossprod(Xmat, y - Z %*% solve(K, crossprod(Z, y))) 
  b <- solve(XVX, XVty) # don't need to scale `XVX / comp[2]` and `XVty / comp[2]

  # compute residuals `r`, r' Vi r & the scaling variance scalar `s2`
  Rmat <- y - Xmat %*% b # r <- as.numeric(Rmat)
  ZtR <- crossprod(Z, Rmat)
  yPy <- crossprod(Rmat) - crossprod(ZtR, solve(K, ZtR)) # 
  s2 <- yPy / comp[2] / nk

  # compute log-likelihood
  ll <- -0.5*nk*(log(2*pi*s2) + 1) - 0.5*Sigma_det_log
  
  if(REML) {
    # need to account for: `XVX` variable = XVX / comp[2] 
    det <- determinant(XVX, log = TRUE)
    log_det_XVX <- as.numeric(det$modulus) - k*log(comp[2])
    
    #log_det_XX <- determinant(crossprod(X), log = TRUE)

    ll <- ll - 0.5*as.numeric(log_det_XVX)
  }

  if(verbose > 1) {
    cat("  -- biglr_ll: gamma", gamma, "; ll", ll, "\n")
  }

  return(as.numeric(ll))
}

biglr_ll_modular <- function(gamma, y, Xmat, Z, K = NULL, REML = TRUE, verbose = 0)
{
  if(verbose > 1) {
    cat("  -- biglr_ll: ")
  }

  if(!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  
  n <- length(y)
  k <- ncol(Xmat)

  nk <- `if`(REML, n - k, n)

  comp <- c(gamma, 1 - gamma)

  # pre-compute `K = Z'Z`
  if(is.null(K)) {
    K <- crossprod(Z)
  }
  
  # compute `Sigma_det_log` = det(V), where V = comp[1] ZZ' + comp[2] I
  Sigma_det_log <- biglr_det(comp, Z, K, log = TRUE)
  
  # compute `XVX` & coeff. `b`, where `XVX` = X' Vi X
  XVt <- biglr_cprodMatInv(comp, Z, Xmat, K, transpose = TRUE) # Vi' X
  XVX <- crossprod(XVt, Xmat) 

  b <- solve(XVX, crossprod(XVt, y)) # (XVX)_inv XV y

  # compute residuals `r`, r' Vi r & the scaling variance scalar `s2`
  Rmat <- y - Xmat %*% b # r <- as.numeric(Rmat)
  rVt <- biglr_cprodMatInv(comp, Z, Rmat, K, transpose = TRUE) # Vi' r
  yPy <- crossprod(rVt, Rmat) 

  s2 <- yPy / nk

  # compute log-likelihood
  ll <- -0.5*nk*(log(2*pi*s2) + 1) - 0.5*Sigma_det_log
  if(REML) {
    log_det_XVX <- determinant(XVX, log = TRUE)
    #log_det_XX <- determinant(crossprod(X), log = TRUE)

    ll <- ll - 0.5*as.numeric(log_det_XVX$modulus)
  }

  if(verbose > 1) {
    cat("gamma", gamma, "; ll", ll, "\n")
  }

  return(as.numeric(ll))
}

naive_ll <- function(gamma, y, X, Z, K = NULL, REML = TRUE)
{
  n <- length(y)
  k <- ncol(X)

  nk <- ifelse(REML, n - k, n)

  G <- tcrossprod(Z)
  Sigma <- gamma*G + (1 - gamma)*diag(n)
  Sigma_inv <- solve(Sigma)
  Sigma_det_log <- as.numeric(determinant(Sigma, log = TRUE)$modulus)

  XV <- crossprod(X, Sigma_inv)
  XVX <- XV %*% X
  b <- solve(XVX) %*% (XV %*% y)

  r <- as.numeric(y - X %*% b)
  yPy <- crossprod(r, Sigma_inv) %*% r
  s2 <- yPy / nk

  ll <- -0.5*nk*(log(2*pi*s2) + 1) - 0.5*Sigma_det_log

  if(REML) {
    log_det_XVX <- determinant(XVX, log = TRUE)
    #log_det_XX <- determinant(crossprod(X), log = TRUE)
  
    ll <- ll - 0.5*as.numeric(log_det_XVX$modulus)
  }

  return(as.numeric(ll))
}
