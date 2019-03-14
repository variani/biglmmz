context("biglmm")

test_that("biglmm: recover true h2", {
  N <- 1000; M <- 200; h2 <- 0.8
  
  Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5
  
  # scale `Zg` for `y` simulation
  col_means <- colMeans(Zg, na.rm = TRUE)
  col_freq <- col_means / 2  # col_means = 2 * col_freq
  col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  
  Z <- sweep(Zg, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")
  
  b <- rnorm(M, 0, sqrt(h2/M))
  y <- Z %*% b + rnorm(N, 0, sqrt(1 - h2))
  
  mod <- biglmm(y, Z = Zg, scale = TRUE)
  
  expect_true(mod$gamma > 0.7)
})

test_that("biglmm: scale Z", {
  N <- 500; M <- 10; h2 <- 0.8
  
  Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5

  # Z is scaled such a way that ZZ' = GRM
  col_means <- colMeans(Zg, na.rm = TRUE)
  col_freq <- col_means / 2  # col_means = 2 * col_freq
  col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  
  Z <- sweep(Zg, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")

  b <- rnorm(M, 0, sqrt(h2/M))
  y <- Zg %*% b + rnorm(N, 0, sqrt(1 - h2))
  
  Zgrm <- Z / sqrt(M) 
  mod_unsc <- biglmm(y, Z = Zgrm)
  mod_sc <- biglmm(y, Z = Zg, scale = TRUE)
  
  expect_equal(mod_unsc$gamma, mod_sc$gamma, tolerance = 1e-8)
})

