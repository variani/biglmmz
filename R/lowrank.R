
#-----------------------
# det = X'Vi
#-----------------------

# Computes det(V)
#' - V = comp[1]*G + comp[2]*diag(N)
#' - G = tcrossprod(Z)
#' - G is a low-rank matrix: ncol(Z) << nrow(Z)
#'
#' (1) \url{https://en.wikipedia.org/wiki/Matrix_determinant_lemma}
#' |A + UCV| = |L| |C| |A|
#' where L = C- + VA-U
#' (2) For our specific case:
#' |D + ZHZ'| = |L| |H| |D|
#' D = diag(comp1), H = diag(comp2), L = (H- + Z'D-Z)
biglr_det <- function(comp, Z, K = NULL, log = TRUE)
{
  stopifnot(log)

  n <- nrow(Z)
  k <- ncol(Z)

  #L <- diag(k) / comp[1] + big_crossprodSelf(Z)[] / comp[2]
  #log_det_L <- determinant(L, log = TRUE)
  #as.numeric(log_det_L$modulus) + n*log(comp[2]) + k*log(comp[1])

  # L = K / comp[2]  or K = comp[2] L
  # det(L) = det(K) / (comp[2])^k
  # logdet(L) = logdet(K) - k log(comp[2])
  if(is.null(K)) {
    K <- crossprod(Z)
  }
  diag(K) <- diag(K) + comp[2] / comp[1]
  log_det_K <- determinant(K, log = TRUE)
  as.numeric(log_det_K$modulus) + n*log(comp[2]) + k*log(comp[1]) - k*log(comp[2])
}


#-----------------------
# cprodMatInv = X'Vi
#-----------------------

#' Computes X'Vi
#'
#' Computes X'Vi in a optimized way given that
#'  - Vi is the inverse of matrix V
#'  - V = comp[1]*G + comp[2]*diag(N)
#'  - G = tcrossprod(Z)
#'  - G is a low-rank matrix: ncol(Z) << nrow(Z)
#'
#' @param comp two-element vector of variance components
#' @param Z FBM of scaled genotypes
#' @param X matrix of covariates
#' @param transpose (default FALSE) logical
#'        whether return transposed product (X'Vi)' = Vi' X
#'
#' @return a product X'Vi
biglr_cprodMatInv <- function(comp, Z, Xmat, K = NULL, transpose = FALSE)
{
  # n <- nrow(Z)
  # k <- ncol(Z)

  # (1) Li <- solve(diag(k) / comp[1] + crossprod(Z) / comp[2])
  # (2) t(X) / comp[2] - crossprod(X, Z) %*% tcrossprod(Li, Z) / (comp[2] * comp[2])

  # L = K / comp[2]  or K = comp[2] L
  if(is.null(K)) {
    K <- crossprod(Z)
  }
  diag(K) <- diag(K) + comp[2] / comp[1]

  part <- Z %*% solve(K, crossprod(Z, Xmat))

  if (transpose) {
    (Xmat - part) / comp[2]
  } else {
    t(Xmat - part) / comp[2]
  }
}

#' Computes X'Vi in an efficient way in a low-rank scenario (see below)
#'
#'  - Vi is the inverse of matrix V
#'  - V = comp[1]*G + comp[2]*diag(N)
#'  - G = tcrossprod(Z)
#' A low-rank scenario: G is a low-rank matrix,
#' i.e. ncol(Z) << nrow(Z)
lr_cprodMatInv <- function(comp, Z, X)
{
  n <- nrow(Z)
  k <- ncol(Z)

  Li <- solve(diag(k) / comp[1] + crossprod(Z) / comp[2])
  t(X) / comp[2] - crossprod(X, Z) %*% tcrossprod(Li, Z) / (comp[2] * comp[2])
}

#' Computes X'Vi in a naive way
#'
#'  - Vi is the inverse of matrix V
#'  - V = comp[1]*G + comp[2]*diag(N)
#'  - G = tcrossprod(Z)
#' The focal scenario: G is a low-rank matrix,
#' i.e. ncol(Z) << nrow(Z)
naive_cprodMatInv <- function(comp, Z, X)
{
  n <- nrow(Z)
  crossprod(X, solve(comp[1]*tcrossprod(Z) + comp[2]*diag(n)))
}

lr_cprodMatInv2 <- function(comp, Z, X)
{
  n <- nrow(Z)
  k <- ncol(Z)

  K <- big_crossprodSelf(Z)
  K[cbind(1:k, 1:k)] <- K[cbind(1:k, 1:k)] + comp[1] / comp[2]

  part <- big_prodMat(Z, solve(K[], big_cprodMat(Z, X)))

  t(X - part) / comp[2]
}

biglr_cprodMatInv_inv <- function(comp, Z, Xmat)
{
  n <- nrow(Z)
  k <- ncol(Z)
  # p <- ncol(Xmat)

  # Step 1
  #Li <- solve(diag(k) / comp[1] + crossprod(Z) / comp[2])
  Li <- big_crossprodSelf(Z)

  big_apply(Li, function(X, ind, comp) {
    X[, ind] <- diag(length(ind)) / comp[1] + X[, ind] / comp[2]
    NULL
  }, a.combine = 'c', comp = comp)
  Li <- solve(Li[]) # here we convert to a matrix by `Li[]`

  # Step 2
  #t(X) / comp[2] - crossprod(X, Z) %*% tcrossprod(Li, Z) / (comp[2] * comp[2])
  #X' / comp[2] - X'Z Li Z' / (comp[2] * comp[2])
  #X' / comp[2] - X'Z Li Z' / (comp[2] * comp[2])
  #X' / comp[2] - [(Z'X)'] Li Z' / (comp[2] * comp[2])
  # (1) P1 = Z'X # p x k
  #X' / comp[2] - [P1' Li] Z' / (comp[2] * comp[2])
  # (2) P2 = P1' Li # p x k
  #X' / comp[2] - P2 Z' / (comp[2] * comp[2])
  # (3) P2 <- P2[] # p x k
  #X' / comp[2] - [Z P2']' / (comp[2] * comp[2])
  # (4) P3 = Z P2' # n x p
  #X' / comp[2] - P3' / (comp[2] * comp[2])
  # (5) P3 <- P3[] # p x n
  #[X / comp[2] - P3 / (comp[2] * comp[2])]'
  # (6) Pt = X / comp[2] - P3 / (comp[2] * comp[2])
  # (7) P = Pt'
  P1mat <- big_cprodMat(Z, Xmat) # p x k
  P2mat <- crossprod(P1mat, Li) # p x k
  P3mat <- big_prodMat(Z, t(P2mat)) # p x n

  # Pt
  Pt <- Xmat / comp[2] - P3mat / (comp[2] * comp[2])
  t(Pt)
}

biglr_cprodMatInv_inv_zt <- function(comp, Zt, Xmat)
{
  n <- ncol(Zt)
  k <- nrow(Zt)
  # p <- ncol(Xmat)

  # Step 1
  #Li <- solve(diag(k) / comp[1] + tcrossprod(Zt) / comp[2])
  Li <- big_tcrossprodSelf(Zt)

  big_apply(Li, function(X, ind, comp) {
    X[, ind] <- diag(length(ind)) / comp[1] + X[, ind] / comp[2]
    NULL
  }, a.combine = 'c', comp = comp)
  Li <- solve(Li[]) # here we convert to a matrix by `Li[]`

  # Step 2
  #X' / comp[2] - X' Zt' Li Zt / (comp[2] * comp[2])
  #X' / comp[2] - (Zt X)' Li Zt / (comp[2] * comp[2])
  #X' / comp[2] - P1' Li Zt / (comp[2] * comp[2])
  # (1) P1 = Zt X # k x p
  #X' / comp[2] - P2  Zt / (comp[2] * comp[2])
  # (2) P2 = P1' Li # p x k
  # Pt = X / comp[2] - Zt' P2' / (comp[2] * comp[2])
  # (3) P = Pt'
  P1mat <- big_prodMat(Zt, X)
  P2mat <- crossprod(P1mat, Li)
  P3mat <- big_cprodMat(Zt, t(P2mat))

  Pt <- X / comp[2] - P3mat / (comp[2] * comp[2])
  t(Pt)
}

biglr_cprodMatInv_evd <- function(comp, Z, Xmat)
{
  n <- nrow(Z)
  k <- ncol(Z)
  # p <- ncol(Xmat)

  # the product
  # L <- diag(k) / comp[1] + crossprod(Z) / comp[2]
  # P = t(X) / comp[2] - X' Z Li Z' / (comp[2] * comp[2])
  #
  # (SVD) Z = U S V'
  # (EVD) Z'Z = V D V', where D = S^2
  # (EVD) L = V K V' = V [(1/comp[1]) I + (1/comp[2]) diag(s^2)] V'
  # (inv) Li = V K^{-1} V'
  #(prod) Z Li Z' = Z V K^{-1} V' Z = W W',
  #       where W = Z V K^{-0.5}
  # Finally:
  # X' Z Li Z' = X' W W' = [X' W] W' = [W' X]' W' = (W [W' X])' = (W [W' X]')'
  #
  # P = [X / comp[2] - W W' X / (comp[2] * comp[2])]'

  # doesn't work: svd <- big_SVD(Z, k = ncol(Z))
  # but this works: svd <- big_SVD(Z, k = ncol(Z) - 1)
  ZZ <- big_crossprodSelf(Z)
  evd <- eigen(ZZ[])

  V <- evd$vectors
  d <- evd$values

  k <- (1/comp[1]) + (1/comp[2]) * d
  W <- big_prodMat(Z, V %*% diag(1/sqrt(k)))

  # v1
  #t(X) / comp[2] - X' W W' / (comp[2] * comp[2])
  t(X) / comp[2] - crossprod(X, W) %*% t(W) / (comp[2] * comp[2])

  # v2
  #t(X / comp[2] - tcrossprod(W, crossprod(W, X))  / (comp[2] * comp[2]))
}


