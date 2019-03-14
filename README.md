
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biglmmz

[![travis-ci build
status](https://travis-ci.org/variani/biglmmz.svg?branch=master)](https://travis-ci.org/variani/biglmmz)

Low-rank linear mixed models (LMMs) powered by
[bigstatsr](https://github.com/privefl/bigstatsr).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("variani/biglmmz")
```

## Example

This is a sanity check to see whether `biglmmz` recovers the
heritability of a quantitative trait on simulated data (1500
individuals, 200 genetic markers, 80% heritability).

``` r
library(biglmmz)

N <- 1500; M <- 200; h2 <- 0.8

# simulate genotypes
set.seed(33)
freqs <- rep(0.5, M) # allele freq. = 0.5
Z <- sapply(freqs, function(f) rbinom(N, 2, f)) 

# scale genotypes
Z_means <- colMeans(Z, na.rm = TRUE)
Z_freq <- Z_means / 2  # Z_means = 2 * Z_freq
Z_sd <- sqrt(2 * Z_freq * (1 - Z_freq))

Z_sc <- sweep(Z, 2, Z_means, "-")
Z_sc <- sweep(Z_sc, 2, Z_sd , "/")

b <- rnorm(M, 0, sqrt(h2/M))
y <- Z_sc %*% b + rnorm(N, 0, sqrt(1 - h2))

# fit model on raw genotypes
m1 <- biglmmz(y, Z = Z, scale = TRUE)
m1$gamma
#> [1] 0.7963861

# fit model on scaled genotypes and normalized by sqrt(M)
Z_norm <- Z_sc / sqrt(M)
m2 <- biglmmz(y, Z = Z_norm, scale = FALSE)
m2$gamma
#> [1] 0.7963861
```
