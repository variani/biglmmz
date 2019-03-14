context("biglr_ll")

test_that("biglr_ll", {
  N <- 5e2; M <- 1e2; P <- 10
  
  y <- rnorm(N)
  Xmat <- matrix(rnorm(P*N), nrow = N, ncol = P)
  comp <- c(0.8, 0.2)
  gamma <- comp[1] / sum(comp)
  
  Z <- FBM(N, M, init = rnorm(N*M))
  Zmat <- Z[]
  
  # REML
  inline_reml <- biglr_ll_inline(gamma, y, Xmat, Z, REML = TRUE)
  modular_reml <- biglr_ll_modular(gamma, y, Xmat, Z, REML = TRUE)
  naive_reml <- naive_ll(gamma, y, Xmat, Zmat, REML = TRUE)

  expect_equal(naive_reml, inline_reml, tolerance = 1e-8)
  expect_equal(naive_reml, modular_reml, tolerance = 1e-8)

  # ML
  modular_ml <- biglr_ll_modular(gamma, y, Xmat, Z, REML = FALSE)
  inline_ml <- biglr_ll_inline(gamma, y, Xmat, Z, REML = FALSE)
  naive_ml <- naive_ll(gamma, y, Xmat, Zmat, REML = FALSE)

  expect_equal(naive_ml, inline_ml, tolerance = 1e-8)  
  expect_equal(naive_ml, modular_ml, tolerance = 1e-8)
})

test_that("biglr_ll: internal", {
  N <- 5e2; M <- 1e2; P <- 10
  
  Xmat <- matrix(rnorm(P*N), nrow = N, ncol = P)
  comp <- c(0.8, 0.2)
  
  Z <- FBM(N, M, init = rnorm(N*M))
  Zmat <- Z[]
  
  Vmat <- comp[1]*tcrossprod(Zmat) + comp[2]*diag(N)
  Vmat_inv <- solve(Vmat)
  
  # log(det(V))
  det_naive <- determinant(Vmat, log = TRUE)
  logdet_naive <- as.numeric(det_naive$modulus)
  
  logdet_biglr <- biglr_det(comp, Z, log = TRUE)
  
  expect_equal(logdet_naive, logdet_biglr, tolerance = 1e-8)
 
  # X'VX
  XVX_naive <- crossprod(Xmat, Vmat_inv) %*% Xmat
  
  XVt <- biglr_cprodMatInv(comp, Z, Xmat, transpose = TRUE) 
  XVX_biglr <- crossprod(XVt, Xmat)
   
  expect_true(all.equal(XVX_naive, XVX_biglr, tolerance = 1e-8))
})

