context("cprodMatIn")

test_that("biglr_cprodMatInv", {
  N <- 5e2; M <- 1e2; P <- 10
  
  Xmat <- matrix(rnorm(P*N), nrow = N, ncol = P)
  comp <- c(0.8, 0.2)
  
  Z <- FBM(N, M, init = rnorm(N*M))
  Zmat <- Z[]
  
  prod_naive <- naive_cprodMatInv(comp, Zmat, Xmat)
  prod_biglr <- biglr_cprodMatInv(comp, Z, Xmat)
  
  expect_true(all.equal(prod_naive, prod_biglr, tolerance = 1e-8))
})

