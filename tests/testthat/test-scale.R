context("scale")

test_that("scale_Z (FBM); M = ncol(Z)", {
  # X => scaled(X) such that X'X = Xsc Xsc / M = GRM
  X <- FBM(5, 5, init = 1:25)
  Xsc1 <- scale(X[]) / sqrt(ncol(X))

  scale_Z(X, M = ncol(X))
  Xsc2 <- X[]

  expect_equal(as.numeric(Xsc1), as.numeric(Xsc2))
})

test_that("scale_Z (FBM); M = 1", {
  X <- FBM(5, 5, init = 1:25)
  Xsc1 <- scale(X[])

  scale_Z(X, M = 1)
  Xsc2 <- X[]

  expect_equal(as.numeric(Xsc1), as.numeric(Xsc2))
})
