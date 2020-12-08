#' Compute the effective sample size for a LMM.
#'
#' @param G A FBM matrix of genotypes. Missing values are not handled.
#' @param cols A vector of columns of G to be used in the model.
#'   By default, all columns of G are used.
#' @param M A scalar for normalization of the
#'   genetic relationship matrix: GRM = Z'Z / M,
#'   where Z is a matrix of standardized genotypes.
#'  By default, M = length(cols).
#' @param h2 The estimated heritability by LMM.
#' @param s2 The estimated scaling constant for the variance components in LMM.
#' @return A data.frame of results.
#'
#' @examples
#' G <- attach_example200()
#' G
#' G[1:5, 1:10]
#'
#' ess(G, h2 = 0.5)
#'
#' ess(G, h2 = 0.8)
#'
#' @export
ess <- function(G, cols = seq(ncol(G)), M = length(cols), h2, s2 = 1.0)
{ 
  ## check arguments
  stopifnot(!missing(G))
  stopifnot(!missing(h2))

  stopifnot(inherits(G, "FBM"))

  ## compute K
  f_sc <- big_scale_grm(M = M) # M is defined in arguments
  stats <- f_sc(G, ind.col = cols)

  K <- big_crossprodSelf(G, fun.scaling = f_sc, ind.col = cols)[]

  # EVD of K
  lamdas <- eigen(K)$values

  # the multiplier
  #  - varianace components = s2 * c(h2, 1 - h2)
  N <- nrow(G)
  mult <- (1/s2) * (sum(1/(h2*lamdas + (1-h2))) + (N-M)/(1-h2)) / N

  ## return results
  res <- data.frame(N = N, M = M, h2 = h2, s2 = s2, mult = mult,
    ESS = N*mult)
  res
}
