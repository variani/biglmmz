context("biglr_det")

test_that("biglr_det", {
  N <- 5e2; M <- 1e2; P <- 10
  
  comp <- c(0.8, 0.2)
  Z <- FBM(N, M, init = rnorm(N*M))
  Zmat <- Z[]
  
  Vmat <- comp[1]*tcrossprod(Zmat) + comp[2]*diag(N)
   
  det_naive <- determinant(Vmat, log = TRUE)
  logdet_naive <- as.numeric(det_naive$modulus)
  
  logdet_biglr <- biglr_det(comp, Z, log = TRUE)
  
  expect_equal(logdet_naive, logdet_biglr, tolerance = 1e-8)
})
