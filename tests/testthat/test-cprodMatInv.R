context("cprodMatIn")

N <- 50; M <- 10; P <- 4

Xmat <- matrix(rnorm(P*N), nrow = N, ncol = P)
comp <- c(0.8, 0.2)

Z <- FBM(N, M, init = rnorm(N*M))
Zmat <- Z[]
Zmat_sc <- scale(Zmat)
Zmat_grm <- Zmat_sc / sqrt(M)

test_that("biglr_cprodMatInv", {
  prod_naive <- naive_cprodMatInv(comp, Zmat, Xmat)
  prod_biglr <- biglr_cprodMatInv(comp, Z, Xmat)
  
  expect_true(all.equal(prod_naive, prod_biglr, tolerance = 1e-8))
})

test_that("biglr_cprodMatInv2 (M = 1)", {
  prod_naive <- naive_cprodMatInv(comp, Zmat_sc, Xmat)
  prod2_biglr <- biglr_cprodMatInv2(comp, Z, Xmat, M = 1)
  
  expect_true(all.equal(prod_naive, prod2_biglr, tolerance = 1e-8))
})

test_that("biglr_cprodMatInv2 (M = ncol(Z)", {
  prod_naive <- naive_cprodMatInv(comp, Zmat_grm, Xmat)
  prod2_biglr <- biglr_cprodMatInv2(comp, Z, Xmat)
  
  expect_true(all.equal(prod_naive, prod2_biglr, tolerance = 1e-8))
})

