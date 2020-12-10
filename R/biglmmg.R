#------------------------
# Main function biglmmg
#------------------------

#' Fit a low-rank LMM on raw genotypes (no explicit scaling).
#'
#' @param y A vector of trait values (quantitative trait).
#' @param X A matrix of covariates. The default value is matrix of ones
#'   for the intercept (one column).
#' @param G A FBM matrix of genotypes. Missing values are not handled.
#' @param cols A vector of columns of G to be used in the model.
#'   By default, all columns of G are used.
#' @param M A scalar for normalization of the
#'   genetic relationship matrix: GRM = Z'Z / M,
#'   where Z is a matrix of standardized genotypes.
#'  By default, M = length(cols).
#' @param K A matrix with the pre-computed cross-product Z'Z / M.
#'  By default, K = NULL, that means K is pre-computed inside the function.
#' @param REML A boolean specifying the likelihood function, REML or ML.
#' @param compute_mult A boolean enabling the computation of 
#'   the effective sample size, when only model fitting is needed.
#'   The default value is TRUE.
#' @param verbose The verbose level.
#'   The default value is 0 (verbose).
#'
#' @details
#' The linear mixed model (LMM) is:
#' y_i =  X_i b + u_i + e_i, where
#'
#' u ~ N(0, s2*h2*G) and e ~ N(0, s2 I)
#'
#' var(y) = V = s2 * (h2*G + I)
#'
#' @examples
#' G <- attach_example200() # load simulated genotypes
#' G
#' G[1:5, 1:10]
#'
#' y <- rnorm(nrow(G)) # generate a random phenotype
#' head(y)
#'
#' mod <- biglmmg(y, G = G)
#' mod$gamma # estimated h2
#'
#' @export
biglmmg <- function(y, X,
  G, cols = seq(ncol(G)), M = length(cols),
  K = NULL,
  REML = TRUE,
  compute_mult = TRUE,
  verbose = 0)
{
  if(verbose) { cat("biglmmg:\n") }

  ### args
  stopifnot(!missing(y))
  if(!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  if(missing(X)) {
    X <- matrix(1, nrow = length(y), ncol = 1)
  }
  stopifnot(!missing(G))

  stopifnot(nrow(y) == nrow(G))
  stopifnot(nrow(y) == nrow(X))

  ## prepare scaling function & stata
  f_sc <- big_scale_grm(M = M) # M is defined in arguments
  stats <- f_sc(G, ind.col = cols)

  ### pre-compute K = Z'Z, where Z = scaled(G) / sqrt(M)
  if(is.null(K)) {
    if (verbose) cat(" - precompute K = crossprod(Z)\n")
      K <- big_crossprodSelf(G, fun.scaling = f_sc, ind.col = cols)[]
  } else {
    if (verbose) cat(" - passed by argument K = crossprod(Z)\n")
    stopifnot(ncol(K) == length(cols))
    stopifnot(nrow(K) == length(cols))
  }

  ### optimize
  if (verbose) cat(" - optimize\n")
  out <- optimize(biglr_ll_grm, c(0, 1),
    y = y, X = X,
    G = G, cols, f_sc = f_sc, stats = stats,
    K = K, REML = REML, verbose = verbose,
    maximum = TRUE)

  r2 <- out$maximum
  ll <- out$objective
  convergence <- NA

  ### estimates of fixed effects
  if (verbose) cat(" - estimate fixed effects\n")

  est <- biglr_fixef_grm(r2, y, X,
    G, cols = cols,
    f_sc = f_sc, stats = stats,
    K = K, REML = TRUE)
  coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  coef$z <- coef$beta / coef$se

  ### estimates of random effects
  if (verbose) cat(" - estimate random effects\n")

  gamma <- r2
  s2 <- est$s2
  comp <- s2 * c(gamma, 1 - gamma)

  # trace factor
  trace_factor <- mult <- NA
  if(compute_mult) {
    if (verbose) cat(" - multiplier\n")
    lamdas <- eigen(K)$values

    N <- length(y)
    M <- nrow(K)

    trace_factor <- (sum(1/(gamma*lamdas + (1-gamma))) + (N-M)/(1-gamma)) / N
    mult <- (1/s2) * trace_factor
  }

  ### return
  mod <- list(
    gamma = gamma, s2 = s2, comp = comp,
    trace_factor = trace_factor, mult = mult,
    est = est, coef = coef,
    REML = REML)

  mod$lmm <- list(r2 = r2, ll = ll, convergence = convergence, REML = REML)

  return(mod)
}

#--------------------------------------------
# Log-likelihood function for low-rank LMM
#--------------------------------------------

biglr_ll_grm <- function(gamma, y, Xmat,
  G, cols = seq(ncol(G)), M = length(cols),
  f_sc = NULL, stats = NULL,
  K, REML = TRUE, verbose = 0)
{
  stopifnot(!missing(K))
  stopifnot(!is.null(K))

  if(is.null(f_sc)) {
    f_sc <- big_scale_grm(M = M) # M is defined in arguments
  }
  if(is.null(stats)) {
    stats <- f_sc(G, ind.col = cols)
  } else {
    stopifnot(nrow(stats) == M)
  }

  if(!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  n <- length(y)
  k <- ncol(Xmat)
  nk <- `if`(REML, n - k, n)

  ng <- length(cols)

  comp <- c(gamma, 1 - gamma)

  # compute `Sigma_det_log` = det(V), where V = comp[1] ZZ' + comp[2] I
  diag(K) <- diag(K) + comp[2] / comp[1]
  log_det_K <- determinant(K, logarithm = TRUE)
  Sigma_det_log <- as.numeric(log_det_K$modulus) + n*log(comp[2]) +
    ng*log(comp[1] / comp[2])

  # compute `XVX` & coeff. `b`, where `XVX` = X' Vi X
  # ZtX <- crossprod(Z, Xmat)
  ZtX <- big_cprodMat(G, Xmat, ind.col = cols, center = stats$center, scale = stats$scale)
  XVX <- (crossprod(Xmat) - crossprod(ZtX, solve(K, ZtX)))

  # XVty <- crossprod(Xmat, y - Z %*% solve(K, crossprod(Z, y)))
  XVty <- crossprod(Xmat, y -
    big_prodMat(G,
      solve(K, big_cprodMat(G, y, ind.col = cols, center = stats$center, scale = stats$scale)),
      ind.col = cols, center = stats$center, scale = stats$scale))
  b <- solve(XVX, XVty) # don't need to scale `XVX / comp[2]` and `XVty / comp[2]

  # compute residuals `r`, r' Vi r & the scaling variance scalar `s2`
  Rmat <- y - Xmat %*% b # r <- as.numeric(Rmat)
  # ZtR <- crossprod(Z, Rmat)
  ZtR <- big_cprodMat(G, Rmat, ind.col = cols, center = stats$center, scale = stats$scale)
  yPy <- crossprod(Rmat) - crossprod(ZtR, solve(K, ZtR)) #
  s2 <- yPy / comp[2] / nk

  # compute log-likelihood
  ll <- -0.5*nk*(log(2*pi*s2) + 1) - 0.5*Sigma_det_log

  if(REML) {
    # need to account for: `XVX` variable = XVX / comp[2]
    det <- determinant(XVX, logarithm = TRUE)
    log_det_XVX <- as.numeric(det$modulus) - k*log(comp[2])

    #log_det_XX <- determinant(crossprod(X), logarithm = TRUE)

    ll <- ll - 0.5*as.numeric(log_det_XVX)
  }

  if(verbose > 1) {
    cat("  -- biglr_ll2: gamma", gamma, "; ll", ll, "\n")
  }

  return(as.numeric(ll))
}

#-------------------------------
# Estimates of fixed effects
#-------------------------------

biglr_fixef_grm <- function(
  gamma, y, Xmat,
  G, cols = seq(ncol(G)), M = length(cols),
  f_sc = NULL, stats = NULL,
  s2, K = NULL, REML = TRUE)
{
  missing_s2 <- missing(s2)

  if(is.null(f_sc)) {
    f_sc <- big_scale_grm(M = M) # M is defined in arguments
  }
  if(is.null(stats)) {
    stats <- f_sc(G, ind.col = cols)
  } else {
    stopifnot(nrow(stats) == M)
  }

  n <- length(y)
  k <- ncol(Xmat)
  nk <- `if`(REML, n - k, n)

  # pre-compute `K = Z'Z`
  if(is.null(K)) {
    K <- big_crossprodSelf(G, fun.scaling = f_sc, ind.col = cols)[]
  } else {
    stopifnot(ncol(K) == length(cols))
    stopifnot(nrow(K) == length(cols))
  }

  # need to compute `s2`
  # then we can get beta/se taking into account `s2`: V = s2 (gamma ZZ' + (1-gamma) I)
  if(missing_s2) {
    comp <- c(gamma, 1 - gamma)

    # compute `XVX` & coeff. `b`, where `XVX` = X' Vi X
    # XVt <- biglr_cprodMatInv(comp, Z, Xmat, K, transpose = TRUE) # Vi' X
    XVt <- biglr_cprodMatInv_grm(comp,
      G, cols = cols, f_sc = f_sc, stats = stats,
      X = Xmat, K = K, transpose = TRUE) # Vi' X
    XVX <- crossprod(XVt, Xmat)
    b <- solve(XVX, crossprod(XVt, y)) # (XVX)_inv XV y

    Rmat <- y - Xmat %*% b # r <- as.numeric(Rmat)
    # rVt <- biglr_cprodMatInv(comp, Z, Rmat, K, transpose = TRUE) # Vi' r
    rVt <- biglr_cprodMatInv_grm(comp,
      G, cols = cols, f_sc = f_sc, stats = stats,
      X = Rmat, K = K, transpose = TRUE) # Vi' r
    yPy <- crossprod(rVt, Rmat)
    s2 <- as.numeric(yPy / nk)
  }

  # update `comp` from unscaled to scaled by `s2`
  comp <- s2 * c(gamma, 1 - gamma)

  XVt <- biglr_cprodMatInv_grm(comp,
    G, cols = cols, f_sc = f_sc, stats = stats,
    X = Xmat, K = K, transpose = TRUE) # Vi' X
  XVX <- crossprod(XVt, Xmat)
  b <- solve(XVX, crossprod(XVt, y)) # (XVX)_inv XV y
  bcov <- solve(XVX)

  ### return
  out <- list(s2 = s2, b = as.numeric(b), bcov = bcov)
}

#--------------------
# Scaling function
#--------------------

#' Scale and normalize genotypes for GRM.
#'
#' The function creates a scaling function 
#' to be used for linear algebra operations on FBMs
#' through the fun.scaling argument in downstream functions from bigstatsr.
#'
#' The scaling function works such that the input matrix of raw genotypes
#' is represented as a matrix Z, 
#' which cross product is the genetic relationship matrix (GRM),
#' defined as scale(G) scale(G)' / M
#'
#' @param center Center?
#' @param scale Scale?
#' @param M A normalization value. It is typically equal to sqrt(M).
#' @return A function that returns a data.frame with two columns 
#'   "center" and "scale".
#'
#' @export
big_scale_grm <- function(center = TRUE, scale = TRUE, M)
{
  function(X, ind.row = rows_along(X), ind.col = cols_along(X)) {
    m <- length(ind.col)
    if (center) {
      tmp <- big_colstats(X, ind.row, ind.col)
      means <- tmp$sum/length(ind.row)
      sds <- if (scale)
          sqrt(tmp$var)
      else rep(1, m)
    }
    else {
      means <- rep(0, m)
      sds <- rep(1, m)
    }
    data.frame(center = means, scale = sqrt(M)*sds)
  }
}

#---------------------------
# Linear Algebra functions
#---------------------------

biglr_cprodMatInv_grm <- function(comp,
  G, cols = seq(ncol(G)), M = length(cols),
  f_sc = NULL, stats = NULL,
  X, K = NULL, transpose = FALSE)
{
  if(is.null(f_sc)) {
    f_sc <- big_scale_grm(M = M) # M is defined in arguments
  }
  if(is.null(stats)) {
    stats <- f_sc(G, ind.col = cols)
  } else {
    stopifnot(nrow(stats) == M)
  }

  if(is.null(K)) {
    K <- big_crossprodSelf(G, fun.scaling = f_sc, ind.col = cols)[]
  } else {
    stopifnot(ncol(K) == length(cols))
    stopifnot(nrow(K) == length(cols))
  }
  diag(K) <- diag(K) + comp[2] / comp[1]

  part <- big_prodMat(G,
    solve(K, big_cprodMat(G, X, ind.col = cols,
      center = stats$center, scale = stats$scale)),
    ind.col = cols, center = stats$center, scale = stats$scale)

  if (transpose) {
    (X - part) / comp[2]
  } else {
    t(X - part) / comp[2]
  }
}

