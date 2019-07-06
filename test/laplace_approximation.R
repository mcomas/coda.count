library(coda.base)

x = matrix(c(0, 10), nrow = 1)
mu = c(-1)
sigma = matrix(1) * 10
B = ilr_basis(2)
Binv = t(MASS::ginv(B))

prob = coda.count::dlrnm(x, mu, sigma, B, hermite.order = 500)
M = c_posterior_approximation(x, mu, solve(sigma), Binv)

posterior = function(h_){
  sapply(h_, function(h){
    p = composition(h, basis = B)
    dbinom(x[1], sum(x), p[1]) * dnorm(h, mu, sigma)
  }) / prob
}

h = seq(-20, 10, 0.01)
f1 = posterior(h)
f2 = dnorm(h, M[1,2], sqrt(M[1,1]))
ylim = range(c(f1,f2))
xlim = range(h[pmax(f1,f2) > 0.0001])
plot(h, f1, type = 'l', xlim = xlim, ylim = ylim)
points(h, f2, type='l', col = 'red')

