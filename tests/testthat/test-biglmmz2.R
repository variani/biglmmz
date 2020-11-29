context("biglmmz2")

set.seed(1)
N <- 1000; M <- 200; h2 <- 0.8
Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5
Zsc <- scale(Zg)
Zgrm <- Zsc / sqrt(M)

Z <- as_FBM(Zg)
f_sc <- big_scale_grm(M = M)
stats <- f_sc(Z)

stats0 <- big_scale()(Z)

# scale `Zg` for `y` simulation
b <- rnorm(M, 0, sqrt(h2/M))
# Zb <- Zsc %*% b
Zb <- big_prodVec(Z, b, center = stats0$center, scale = stats0$scale)
y <- Zb + rnorm(N, 0, sqrt(1 - h2))

test_that("biglmmg: recover true h2", {
  mod <- biglmmz(y, Z = Zg, scale = TRUE)
  mod2 <- biglmmg(y, G = Z)

  expect_true(mod$gamma > 0.7)
  expect_true(mod2$gamma > 0.7)
})

test_that("biglmmg: fixed effect", {
  mod <- biglmmz(y, Z = Zg, scale = TRUE)
  mod2 <- biglmmg(y, G = Z)

  expect_equal(mod$coef$beta, mod2$coef$beta, tol = 1e-5)
  expect_equal(mod$coef$se, mod2$coef$se, tol = 1e-5)
})
