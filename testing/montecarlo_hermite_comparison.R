library(coda.base)
library(coda.count)
set.seed(1)

mu = c(0,1,0)
sigma = diag(3)
mc.nsim = 1000
B = ilr_basis(4)

X = rlrnm(10, 100, mu, sigma, B)
Z = randtoolbox::sobol(mc.nsim, dim = length(mu), normal = TRUE)

L_mc = c_lrnm_posterior_moments_montecarlo(t(X), B %*% mu, B %*% sigma %*% t(B), t(Z))
L_hermite = c_lrnm_posterior_moments_hermite(t(X), B %*% mu, B %*% sigma %*% t(B), 10)
cbind(t(L_mc$clr_E1),
      0, 0,
      t(L_hermite$clr_E1))

fit_mc = c_lrnm_fit_montecarlo(t(X), t(Z), 0.001, 100)
fit_hermite = c_lrnm_fit_hermite(t(X), 10, 0.001, 100)
fit_mc$clr_mu
fit_hermite$clr_mu


fit_mc$clr_sigma
fit_hermite$clr_sigma

fit_mc$em_iter
fit_hermite$em_iter

fit_lrnm(X, method = 'mc')$sigma
fit_lrnm(X, method = 'hermite')$sigma
fit_lrnm(X, method = 'laplace')$sigma
