# biglmmz

[![travis-ci build status](https://travis-ci.org/variani/biglmmz.svg?branch=master)](https://travis-ci.org/variani/biglmmz)

Low-rank linear mixed models (LMMs) powered by [bigstatsr](https://github.com/privefl/bigstatsr).

## Example

```
library(biglmmz)

N <- 1500; M <- 200; h2 <- 0.8
  
Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5

col_means <- colMeans(Zg, na.rm = TRUE)
col_freq <- col_means / 2  # col_means = 2 * col_freq
col_sd <- sqrt(2 * col_freq * (1 - col_freq))

Z <- sweep(Zg, 2, col_means, "-")
Z <- sweep(Z, 2, col_sd , "/")

b <- rnorm(M, 0, sqrt(h2/M))
y <- Z %*% b + rnorm(N, 0, sqrt(1 - h2))
  
mod <- biglmmz(y, Z = Zg, scale = TRUE)
mod$gamma
```
