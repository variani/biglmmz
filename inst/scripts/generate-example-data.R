library(bigstatsr)

# simulation parameters
N <- 1500 # sample size
M <- 200 # number of genotypes

# simulate genotypeo33)
freqs <- rep(0.5, M) # allele freq. = 0.5
X <- sapply(freqs, function(f) rbinom(N, 2, f)) 

code <- rep(NA_real_, 256)
code[1:3] <- c(0, 1, 2)
G <- FBM.code256(nrow(X), ncol(X), code, 
  init = X, backingfile = "example200")
G$save()
