#
#                                                         Function_Call
#1 lmm_sc_K<-biglmm(dat$y,Z=dat$Z,scale=TRUE,compute_K=FALSE,verbose=2)
#2    lmm_sc<-biglmm(dat$y,Z=dat$Z,scale=TRUE,compute_K=TRUE,verbose=2)
#3 lmm_K<-biglmm(dat$y,Z=dat$Zsc,scale=FALSE,compute_K=FALSE,verbose=2)
#4    lmm<-biglmm(dat$y,Z=dat$Zsc,scale=FALSE,compute_K=TRUE,verbose=2)
#  Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
#1           16.458                0.1            1655.7
#2           14.759                0.1            1654.9
#3           17.135                0.0             548.9
#4           11.571                0.0             427.0

### inc
library(peakRAM)
library(bigstatsr)

### par
N <- 100e3 # individuals
M <- 500  # markers
h2 <- 0.8

### data simulation
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

  Zsc <- Z / sqrt(M)

  list(y = y, Z = as_FBM(Zg), Zsc = as_FBM(Zsc))
}

dat <- sim_data(N, M, h2)
  
### profiling
prof <- peakRAM(
  lmm_sc_K <- biglmm(dat$y, Z = dat$Z, scale = TRUE, compute_K = FALSE, verbose = 2),
  lmm_sc <- biglmm(dat$y, Z = dat$Z, scale = TRUE, compute_K = TRUE, verbose = 2),
  lmm_K <- biglmm(dat$y, Z = dat$Zsc, scale = FALSE, compute_K = FALSE, verbose = 2),
  lmm <- biglmm(dat$y, Z = dat$Zsc, scale = FALSE, compute_K = TRUE, verbose = 2)
)

print(prof)
