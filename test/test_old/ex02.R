if(!exists('PATTERN')) PATTERN = "K_0002-SEED_0001"

PARS = regmatches(PATTERN, regexec("K_(.*?)-SEED_(.*)", PATTERN))[[1]]

K = as.numeric(PARS[2])
SEED = as.numeric(PARS[3])

library(dcount)
library(randtoolbox)

x = 100+1:K
mu = rep(0,K-1)
sigma = diag(1,K-1)

#(P <- dnm(x, mu = mu, sigma = sigma, 100))
#(mu_exp <- m1_dnm(x, mu = mu, sigma = sigma, 100)/P)
#(var_exp <- m2_dnm(x, mu = mu, sigma = sigma, 100)/P)

r0 = expected00(x = x, mu_ilr = mu, sigma_ilr = sigma, order = 1000)

set.seed(SEED)
NSIM = 1000

mu_exp = c(r0[[2]])

Z1 = rnorm(NSIM*length(mu))
Z1 = matrix(Z1, ncol=length(mu))
r1 = expected01(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z1, mu_exp = mu_exp)

Z2 = rnorm(NSIM*length(mu)/2)
Z2 = matrix(Z2, ncol=length(mu))
Z2 = rbind(Z2,-Z2)
r2 = expected01(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z2, mu_exp = mu_exp)

Z3 = sobol(NSIM, dim = length(mu), scrambling = 2, normal = TRUE, seed = SEED)
Z3 = matrix(Z3, ncol=length(mu))
r3 = expected01(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z3, mu_exp = mu_exp)

Z4 = sobol(NSIM/2, dim = length(mu), scrambling = 2, normal = TRUE, seed = SEED)
Z4 = matrix(Z4, ncol=length(mu))
Z4 = rbind(Z4,-Z4)
r4 = expected01(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z4, mu_exp = mu_exp)

r5 = expected_metropolis(x = x, mu_ilr = mu, sigma_ilr = sigma, mu_exp = mu_exp, nsim = NSIM)

r6 = expected_hamiltonian(x = x, mu_ilr = mu, sigma_ilr = sigma, nsim = NSIM)

df = data.frame(
  'var' = c('m1', 'm2'),
  'err.r1' = sapply(2:3, function(i) sum((r0[[i]]-r1[[i]])^2)),
  'err.r2' = sapply(2:3, function(i) sum((r0[[i]]-r2[[i]])^2)),
  'err.r3' = sapply(2:3, function(i) sum((r0[[i]]-r3[[i]])^2)),
  'err.r4' = sapply(2:3, function(i) sum((r0[[i]]-r4[[i]])^2)),
  'err.r5' = sapply(2:3, function(i) sum((r0[[i]]-r5[[i]])^2)),
  'err.r6' = sapply(2:3, function(i) sum((r0[[i]]-r6[[i]])^2)))

save(list = ls(), file = sprintf("ex02/data-K_%04d-SEED_%04d.RData", K, SEED))
