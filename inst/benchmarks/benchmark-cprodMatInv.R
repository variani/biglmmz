# Results of benchmarking for N = 10e3, M = 500
#                   Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
#naive_cprodMatInv         4.049                  0             763.2
#lr_cprodMatInv            0.075                  0              28.7
#biglr_cprodMatInv_inv3    0.189                  0              74.1
#biglr_cprodMatInv_evd4    0.189                  0             104.2

### inc
library(peakRAM)
library(bigstatsr)

### par
N <- 70e3 # individuals
M <- 500  # markers

### variables
#X <- matrix(1, nrow = N, ncol = 1)
X <- matrix(rnorm(10*N), nrow = N, ncol = 10)
comp <- c(0.5, 0.5)

### simulate Z
Z <- FBM(N, M, init = rnorm(N*M))

Zt <- big_transpose(Z)

Zmat <- Z[]

### profiling
prof <- peakRAM(
  #naive_cprodMatInv(comp, Zmat, X), # hardly possible for N > 10e3
  lr <- lr_cprodMatInv(comp, Zmat, X),
  lr2 <- lr_cprodMatInv2(comp, Z, X),
  biglr_inv <- biglr_cprodMatInv_inv(comp, Z, X),
  biglr_evd <- biglr_cprodMatInv_evd(comp, Z, X),
  biglr_inv_zt <- biglr_cprodMatInv_inv_zt(comp, Zt, X)
)

stopifnot(all.equal(lr, lr2,          tolerance = 1e-8))
stopifnot(all.equal(lr, biglr_inv,    tolerance = 1e-8))
stopifnot(all.equal(lr, biglr_evd,    tolerance = 1e-8))
stopifnot(all.equal(lr, biglr_inv_zt, tolerance = 1e-8))

print(prof)

### Why time(biglr_cprodMatInv_inv) < time(biglr_cprodMatInv_evd)
# which operation `eigen` or `solve` is more expensive?
ZZ <- big_crossprodSelf(Z)
ZtZt <- big_tcrossprodSelf(Zt)
prof2 <- peakRAM(
  eigen(ZZ[]),
  solve(ZZ[]))

print(prof2)

# Zt vs Z
prof3 <- peakRAM(
  big_crossprodSelf(Z),
  big_tcrossprodSelf(Zt))

print(prof3)
