library(coda.dist)
library(rstan)

m = stan_model('test/nm_stan.stan')

expected_hamiltonian = function(x, mu_ilr, sigma_ilr, nsim){
  K = as.integer(length(x))
  x = matrix(x, ncol = 1)
  mu = matrix(mu_ilr, ncol = 1)
  sigma = sigma_ilr
  B = ilr_basis(K)
  MCONSTANT = lpmultinomial_const(x)
  data = list(K = K, x = x, mu = mu, sigma = sigma, B = B, MCONSTANT = MCONSTANT)
  res = sampling(m, iter = nsim + 1000, chain = 1,
                 warmup = 1000, check_data = FALSE, verbose = FALSE, show_messages = FALSE,
                 data = data)
  rh = matrix(rstan::extract(res)$h, ncol = K-1)
  list(
    0,
    colMeans(rh),
    (t(rh) %*% rh) / NROW(rh) )
}

save(m, expected_hamiltonian, file='test/nm_stan_model.RData')

if(FALSE){
  library(coda.dist)
  library(rstan)
  load('test/nm_stan_model.RData')
  K = 2
  x = 1:K
  mu = rep(0,K-1)
  sigma = diag(1,K-1)
  NSIM = 1000
  expected_hamiltonian(x = x, mu_ilr = mu, sigma_ilr = sigma, nsim = NSIM)
}
