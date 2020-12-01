
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

# fit model on raw genotypes & avoid explicit scaling
G <- as_FBM(Z) # input matrix of genotypes is FBM
m3 <- biglmmg(y, G = as_FBM(Z))
#> function(X, ind.row = rows_along(X), ind.col = cols_along(X)) {
#>     m <- length(ind.col)
#>     if (center) {
#>       tmp <- big_colstats(X, ind.row, ind.col)
#>       means <- tmp$sum/length(ind.row)
#>       sds <- if (scale)
#>           sqrt(tmp$var)
#>       else rep(1, m)
#>     }
#>     else {
#>       means <- rep(0, m)
#>       sds <- rep(1, m)
#>     }
#>     data.frame(center = means, scale = sqrt(M)*sds)
#>   }
#> <environment: 0x7fdb26f4c7d8>
m3$gamma 
#> [1] 0.7974539

# Effective sample size (ESS) multipier
K <- big_crossprodSelf(G, fun.scaling = big_scale_grm(M = M))[]
# K <- crossprod(Z_norm)

lamdas <- eigen(K)$values

# varianace components = s2 * c(h2, 1 - h2)
h2 <- m3$gamma
s2 <- m3$s2

mult <- (1/s2) * (sum(1/(h2*lamdas + (1-h2))) + (N-M)/(1-h2)) / N

res <- data.frame(N = N, M = M, h2_hat = h2, s2 = s2, mult = mult)
res
#>      N   M    h2_hat        s2     mult
#> 1 1500 200 0.7974539 0.9743175 4.416904
```
