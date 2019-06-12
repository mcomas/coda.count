SEED = 1

library(coda.dist)

nodes = 1500
x = c(0,100,100)
mu = rep(0,2)
sigma = matrix(c(1,-0.9,
                 -0.9, 1), nrow=2)

B = ilr_basis(length(x))
mu_max = mvf_maximum(x, mu, sigma, B, 0.00001, max_iter = 1000, prop = 0.8)

f0 = function() expected_hermite(x = x, mu_ilr = mu, sigma_ilr = sigma, order = nodes)
r0 = f0()


set.seed(SEED)
NSIM = 10000


f1 = function(){
  Z = rnorm(NSIM*length(mu))
  Z = matrix(Z, ncol=length(mu))
  expected_montecarlo_04(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z, m1 = r0[[2]], m2 = r0[[3]])
}
r1 = f1()

f2 = function(){
  Z = rnorm( (NSIM/2)*length(mu))
  Z = matrix(Z, ncol=length(mu))
  Z = rbind(Z,-Z)
  expected_montecarlo_04(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z, m1 = r0[[2]], m2 = r0[[3]])
}
r2 = f2()

f3 = function(){
  Z = sobol(NSIM, dim = length(mu), scrambling = 2, normal = TRUE, seed = SEED)
  Z = matrix(Z, ncol=length(mu))
  expected_montecarlo_04(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z, m1 = r0[[2]], m2 = r0[[3]])
}
r3 = f3()

f4 = function(){
  Z = sobol((NSIM/2), dim = length(mu), scrambling = 2, normal = TRUE, seed = SEED)
  Z = matrix(Z, ncol=length(mu))
  Z = rbind(Z,-Z)
  expected_montecarlo_04(x = x, mu_ilr = mu, sigma_ilr = sigma, Z = Z, m1 = r0[[2]], m2 = r0[[3]])
}
r4 = f4()

r0[[2]]-r1[[2]]
r0[[2]]-r2[[2]]
r0[[2]]-r3[[2]]
r0[[2]]-r4[[2]]

# MCMC using metropolis algorithm
f5 = function(){
  expected_metropolis(x = x, mu_ilr = mu, sigma_ilr = sigma, mu_exp = mu_max, nsim = NSIM)
}
r5 = f5()

library(rstan)
load('stan_model.RData')
# MCMC using hamiltonian algorithm
f6 = function(){
  expected_hamiltonian(x = x, mu_ilr = mu, sigma_ilr = sigma, nsim = NSIM)
}
r6 = f6()
