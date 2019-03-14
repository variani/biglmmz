
#-----------------------
# Scaling functions
#-----------------------

impute_mean <- function(Z)
{
  mat_na <- is.na(Z)
  means <- colMeans(Z, na.rm = TRUE)
  
  ind_na <- lapply(seq(ncol(mat_na)), function(col) {
    x <- mat_na[, col]
    rows_na <- which(x)
    if(length(rows_na)) {
      cbind(rows_na, col)
    } else {
      NULL
    }
  }) 
  ind_na <- do.call(rbind, ind_na)
  
  Z[ind_na] <- means[ind_na[, 2]]
  
  Z
}

scale_z <- function(Z, impute = FALSE)
{
  if(impute) {
    Z <- impute(Z)
  }
  
  scale(Z, center = TRUE, scale = TRUE)
}

scale_zg <- function(Z, impute = FALSE)
{
  if(impute) {
    Z <- impute(Z)
  }
  
  scale(Z, center = TRUE, scale = TRUE) / sqrt(ncol(Z))
}



