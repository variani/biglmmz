context("coef")

test_that("biglmmz: coef", {
  if(requireNamespace("lme4qtl", quitely = TRUE)) {
    N <- 200; M <- 50; h2 <- 0.8
  
    Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5
  
    # scale `Zg` for `y` simulation
    col_means <- colMeans(Zg, na.rm = TRUE)
    col_freq <- col_means / 2  # col_means = 2 * col_freq
    col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  
    Z <- sweep(Zg, 2, col_means, "-")
    Z <- sweep(Z, 2, col_sd , "/")
  
    b <- rnorm(M, 0, sqrt(h2/M))
    y <- 5 + Z %*% b + rnorm(N, 0, sqrt(1 - h2))
  
    # prepare data & fit models
    ids <- as.character(seq(N))
    Zgrm <- Z / sqrt(M)
    G <- tcrossprod(Zgrm)
    rownames(G) <- colnames(G) <- ids
  
    dat <- data.frame(y = y, id = ids)
    mod0 <- lme4qtl::relmatLmer(y ~ (1|id), dat, relmat = list(id = G))
  
    mod <- biglmmz(y, Z = Zg, scale = TRUE)
  
    expect_true(all.equal(as.numeric(summary(mod0)$coef), as.numeric(mod$coef), tol = 1e-4))
  }
})


