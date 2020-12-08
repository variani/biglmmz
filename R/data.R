#' Attach example genotypes in the FBM format (bigstatsr).
#'
#' The function attaches simulated genotype matrix
#' with 1,500 samples (rows) and 200 genetic variants (colums).
#' Genotypes are coded as 0, 1 and 2 (the number of alternative alleles).
#' All variants have the minor allele frequency 0.5.
#'
#' @return object of class FBM. 
#' @examples
#' G <- attach_example200()
#' G
#' G[1:5, 1:10]
#'
#' @export
attach_example200 <- function()
{
  rds <- system.file("extdata", "example200.rds", package = "biglmmz")
  big_attach(rds)
}