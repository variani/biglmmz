#' Attach example genotypes in the FBM format (bigstatsr).
#'
#' The matrix of 200 genotypes coded as 0, 1 and 2 (columns)
#' in 1,000 individuals (rows).
#'
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
