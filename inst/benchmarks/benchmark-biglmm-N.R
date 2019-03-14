# For vals_N <- c(10e3, 50e3, 100e3:
#
#                                      Function_Call Elapsed_Time_sec
#1 with(dat[[1]],biglmm(y,Z=Z,scale=TRUE,verbose=2))            5.261
#2 with(dat[[2]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           21.684
#3 with(dat[[3]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           42.152
#  Total_RAM_Used_MiB Peak_RAM_Used_MiB
#1                0.2             663.2
#2                0.0             573.4
#3                0.0            1146.8

### inc
library(peakRAM)
library(bigstatsr)
library(biglmm)

### par
vals_N <- c(10e3, 50e3, 100e3) # 300e3 is too extreme for my desktop
M <- 500  # markers
h2 <- 0.8

### function to simulate data
sim_data <- function(N, M, h2)
{
  Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5

  col_means <- colMeans(Zg, na.rm = TRUE)
  col_freq <- col_means / 2  # col_means = 2 * col_freq
  col_sd <- sqrt(2 * col_freq * (1 - col_freq))

  Z <- sweep(Zg, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")

  b <- rnorm(M, 0, sqrt(h2/M))
  y <- Z %*% b + rnorm(N, 0, sqrt(1 - h2))

  list(y = y, Z = as_FBM(Zg))
}

dat <- lapply(vals_N, function(N) sim_data(N, M, h2))

### profiling
prof <- peakRAM(
  with(dat[[1]], biglmm(y, Z = Z, scale = TRUE, verbose = 2)),
  with(dat[[2]], biglmm(y, Z = Z, scale = TRUE, verbose = 2)),
  with(dat[[3]], biglmm(y, Z = Z, scale = TRUE, verbose = 2)))

print(prof)

#### Before modif:
#                                       Function_Call Elapsed_Time_sec
# 1 with(dat[[1]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           11.062
# 2 with(dat[[2]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           45.821
# 3 with(dat[[3]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           92.037
#   Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1                0.1             887.6
# 2                0.0             710.0
# 3                0.0            1146.8
#### After modif:
#                                       Function_Call Elapsed_Time_sec
# 1 with(dat[[1]],biglmm(y,Z=Z,scale=TRUE,verbose=2))            1.911
# 2 with(dat[[2]],biglmm(y,Z=Z,scale=TRUE,verbose=2))            9.946
# 3 with(dat[[3]],biglmm(y,Z=Z,scale=TRUE,verbose=2))           20.412
#   Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1                0.1             327.5
# 2                0.0             853.0
# 3                0.0            1677.4
