SEED = 2
library(coda.dist)
library(normalmultinomial)
library(randtoolbox)
library(dplyr)

x = c(1,0)
MU.orig = 0
SIGMA.orig = matrix(1)

set.seed(SEED)
X = rnormalmultinomial(1000, sample(c(10,100,10000), 100, replace = TRUE), MU.orig, SIGMA.orig)

#### Quasi-random normal sample
NSIM = 10000
Z0 = sobol(NSIM/2, dim = length(MU.orig), normal = TRUE, init=TRUE)
Z = matrix(0, nrow=NSIM, ncol=1)
Z[seq(1,NSIM,2),] = Z0
Z[seq(2,NSIM,2),] = -Z0


H.dm = ilr_coordinates(dm_fit(X)$expected)


H.mu = H.dm
MU = colMeans(H.dm)
SIGMA = cov(H.dm)

# First fit
fit0 = sapply(1:nrow(X), function(i){
  E = expected_montecarlo_01(x = X[i,], mu_ilr = MU, sigma_ilr = SIGMA,
                             Z = Z, mu_exp = H.mu[i,])
  c(E[[2]][1],E[[3]][1])
}) %>% t()

####
fit = fit0
iter = 0
# cbind(fit, fit[,2] - fit[,1] * fit[,1], X)
repeat{
  H.mu <- matrix(fit[,1], ncol = 1)
  H.sigma <- matrix(fit[,2] - fit[,1] * fit[,1], ncol = 1)

  print(MU.next <- colMeans(H.mu))
  iter = iter + 1
  if(abs(MU.next - MU) < 10^-6 | 20 < iter){
    break
  }
  MU = MU.next
  print(iter)
  (SIGMA <- matrix(mean(fit[,2]) - MU^2))

  fit = sapply(1:nrow(X), function(i){
    E = expected_montecarlo_02(x = X[i,], mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = H.mu[i,], var_exp = H.sigma[i,])
    c(E[1,],E[2,])
  }) %>% t()
}

