#----------------------------------
# Impute functions for FBM
#----------------------------------
impute_Z <- function(Z, impute = "zero")
{
  stopifnot(inherits(Z, "FBM"))
  stopifnot(impute == "zero")

  big_apply(Z, function(Z, ind) {
    Z0_part <- Z[, ind]
    Z0_part[is.na(Z0_part)] <- 0
    Z[, ind] <- Z0_part
    NULL
  })

  return(invisible())
}

#----------------------------------
# Scaling functions for FBM
#----------------------------------
scale_Z <- function(Z, impute = TRUE, M = ncol(Z))
{
  stopifnot(inherits(Z, "FBM"))

  big_apply(Z, function(Z, ind) {
    # scale Z such a way that ZZ' = GRM if M = ncol(Z)
    Z0_part <- Z[, ind]

    # col_means <- colMeans(Z0_part, na.rm = TRUE)
    col_means <- matrixStats::colMeans2(Z0_part, na.rm = TRUE)
    # col_freq <- col_means / 2  # col_means = 2 * col_freq
    # col_sd <- sqrt(2 * M * col_freq * (1 - col_freq))
    col_sd <- matrixStats::colSds(Z0_part, center = col_means, na.rm = TRUE)

    # Z0_part <- sweep(sweep(Z0_part, 2, col_means, "-"), 2, col_sd , "/")
    Z0_part <- scale(Z0_part, center = col_means, scale = col_sd * sqrt(M))
    
    if(impute) {
      Z0_part[is.na(Z0_part)] <- 0.0
    }

    Z[, ind] <- Z0_part

    NULL
  })

  return(invisible())
}

#----------------------------------
# Scaling functions for matrices
#----------------------------------

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



