library(coda.count)
X = simplex_lattice(3,3)
mu = c(0,1)
sigma = diag(2)
cbind(X, t(apply(X, 1, c_m1_hermite, mu, sigma, eps = 1e-10)))

