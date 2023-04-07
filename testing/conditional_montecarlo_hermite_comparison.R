library(coda.base)
library(coda.count)
set.seed(1)

mu = c(0,1,0)
sigma = diag(3)
mc.nsim = 10000
B = ilr_basis(4)

X = rlrnm(10, 10, mu, sigma, B)
Z = randtoolbox::sobol(mc.nsim, dim = length(mu), normal = TRUE)

tX = t(X)

str(c_lrnm_posterior_moments_hermite(tX, B %*% mu, B %*% sigma %*% t(B), 10))
str(c_cond_lrnm_posterior_moments_hermite(tX, B %*% mu, B %*% sigma %*% t(B), tX > 0, 10))
str(c_cond_lrnm_posterior_moments_montecarlo(tX, B %*% mu, B %*% sigma %*% t(B), tX > 0, t(Z)))
str(c_cond_lrnm_posterior_moments_hermite(tX, B %*% mu, B %*% sigma %*% t(B), tX > 0, 10))

L_mc = c_lrnm_posterior_moments_montecarlo(t(X), B %*% mu, B %*% sigma %*% t(B), t(Z))
L_hermite = c_lrnm_posterior_moments_hermite(t(X), B %*% mu, B %*% sigma %*% t(B), 10)
cbind(t(L_mc$clr_E1),
      0, 0,
      t(L_hermite$clr_E1))

c_cond_lrnm_fit_montecarlo(X = t(X), C = t(X>0), Z = t(Z), 0.01, 100)

C_free = apply(X, 1, function(x) 1:length(x) == which.max(x))

M1 = cbind(
  c_lrnm_fit_montecarlo(X = t(X), Z = t(Z), 0.01, 100)$clr_mu,
  c_lrnm_fit_hermite(X = t(X), order = 5, 0.0001, 100)$clr_mu,
  c_cond_lrnm_fit_montecarlo(X = t(X), C = C_free, Z = t(Z), 0.01, 100)$clr_mu,
  c_cond_lrnm_fit_hermite(X = t(X), C = C_free, order=5, 0.0001, 100)$clr_mu,
  c_low_dim_cond_lrnm_fit_hermite(X = t(X), C = C_free, maxd0 = 5, order = 5, 0.0001, 100)$clr_mu)
M1

# I would expect that when maxd0 is big low_dim_cond_lrnm = cond_lrnm
I1 = c_lrnm_fit_montecarlo(X = t(X), Z = t(Z), 0.01, 1)
I2 = c_cond_lrnm_fit_montecarlo(X = t(X), C = C_free, Z = t(Z), 0.01, 1)

# Rotation for correlated sample is different in c_lrnm and c_cond_lrnm because of basis
(I1$Zk[,1:10,1] -
    (I1$Bs_N_mu[,1] + t(chol(I1$Bs_N_sigma[,,1])) %*% t(Z)[,1:10])) |>
  range()

(I2$Zk[,1:10,1] -
    (I2$B0s_N_mu[,1] + t(chol(I2$B0s_N_sigma[,,1])) %*% t(Z)[,1:10])) |>
  range()

cov2cor(I1$Bs_N_sigma[,,1])
cov2cor(I2$B0s_N_sigma[,,1] %*% I2$B0s[,,1] %*% t(I2$B0s[,,1]))

I1$Bs
I2$B[,,1]

eigen(I1$Bs_N_sigma[,,1])$vectors
eigen(I1$Bs_N_sigma[,,1])$values
eigen(I2$B0s_N_sigma[,,1])$vectors
eigen(I2$B0s_N_sigma[,,1])$values

t(I1$Bs) %*% I1$Zk[,1:10,1]
t(I2$B[,,1]) %*% t(I2$B0s[,,1]) %*% I2$Zk[,1:10,1]

