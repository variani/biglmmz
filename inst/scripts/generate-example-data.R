library(bigstatsr)

# simulation parameters
N <- 1500 # sample size
M <- 200 # number of genotypes

# simulate genotypes
set.seed(33)
freqs <- rep(0.5, M) # allele freq. = 0.5
X <- sapply(freqs, function(f) rbinom(N, 2, f)) 

G <- as_FBM(X, type = "integer", backingfile = "example200")
G$save()
