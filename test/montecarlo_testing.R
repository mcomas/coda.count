X = rbind(c(1,2,3),
          c(3,2,1))
x = X[1,]
mu = c(0,1)
sigma = matrix(c(1,0.8,0.8,2), nrow = 2)

B = ilr_basis(3)

dlrnm(x, mu, sigma, B, method = 'hermite', 10)
Z = matrix(rnorm(2*1000), ncol = 2)
c_d_lrnm_montecarlo(x, mu, sigma, t(MASS::ginv(B)), Z)
library(randtoolbox)
Z = sobol(100000, dim = 2, normal = TRUE)
c_d_lrnm_montecarlo(x, mu, sigma, t(MASS::ginv(B)), Z)

Z = halton(100, dim = 2, normal = TRUE)
c_d_lrnm_montecarlo(x, mu, sigma, t(MASS::ginv(B)), Z)

c_obtain_moments_lrnm_hermite(X, mu, sigma, B, 20)
c_obtain_moments_lrnm_montecarlo(X, mu, sigma, B, Z)
