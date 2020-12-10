
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biglmmz

<!-- badges: start -->

[![R build
status](https://github.com/privefl/biglmmz/workflows/R-CMD-check/badge.svg)](https://github.com/privefl/biglmmz/actions)
<!-- badges: end -->

Low-rank linear mixed models (LMMs) powered by
[bigstatsr](https://github.com/privefl/bigstatsr).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("variani/biglmmz")
```

## Examples

We first check whether the low-rank linear mixed model (LMM) recovers
the heritability of a quantitative trait on simulated data (1500
individuals, 200 genetic markers, heritability 80%).

``` r
library(biglmmz)

## load simulated genotypes
(G <- attach_example200())
```

    #> A Filebacked Big Matrix of type 'code 256' with 1500 rows and 200 columns.

``` r
(N <- nrow(G))
```

    #> [1] 1500

``` r
(M <- ncol(G))
```

    #> [1] 200

``` r
G[1:5, 1:10]
```

    #>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    #> [1,]    2    1    0    2    1    1    1    1    2     1
    #> [2,]    1    2    1    1    2    1    1    2    1     1
    #> [3,]    2    2    1    1    1    1    1    0    2     0
    #> [4,]    0    2    0    2    1    0    0    1    0     2
    #> [5,]    1    2    0    2    2    1    0    0    1     0

``` r
## simulate a phenotype with heritability h2 = 0.8 
h2 <- 0.8
# generate effect sizes under the polygenic model
b <- rnorm(M, 0, sqrt(h2/M))
# pre-compute mean and sd values for each genotype
stats <- big_scale()(G) 
# compute the matrix-vector product, Z b
Zb <- big_prodVec(G, b, center = stats$center, scale = stats$scale)
y <- Zb + rnorm(N, 0, sqrt(1 - h2))

## fit low-rank linear mixed model
mod <- biglmmg(y, G = G)
# check the estimate of h2
mod$gamma 
```

    #> [1] 0.7972606

We next compute the effective sample size (ESS) multiplier for the LMM.
See the
[pre-print](https://www.biorxiv.org/content/10.1101/2019.12.15.877217v2.full).

Calling the [ess](https://variani.github.io/biglmmz/reference/ess.html)
function:

``` r
# varianace components = s2 * c(h2, 1 - h2)
h2 <- mod$gamma
s2 <- mod$s2

ess(G, h2 = h2, s2 = s2)
```

    #>      N   M        h2       s2    mult     ESS
    #> 1 1500 200 0.7972606 1.048888 4.09906 6148.59

Calculating the ESS manually:

``` r
# EVD on K = Z'Z/M, where Z is a matrix of scaled genotypes G
K <- big_crossprodSelf(G, fun.scaling = big_scale_grm(M = M))[]
# EVD of K
lambdas <- eigen(K, only.values = TRUE)$values

# the multiplier
mult <- (1/s2) * (sum(1/(h2*lambdas + (1-h2))) + (N-M)/(1-h2)) / N

## print results
(res <- data.frame(N = N, M = M, h2_hat = h2, s2 = s2, mult = mult))
```

    #>      N   M    h2_hat       s2    mult
    #> 1 1500 200 0.7972606 1.048888 4.09906
