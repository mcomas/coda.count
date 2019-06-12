source('ex02.R')


B = ilr_basis(length(x))
mu_max = mvf_maximum(x, mu, sigma, B, 0.00001, max_iter = 1000, prop = 0.8)
mu_exp = expected00(x = x, mu_ilr = mu, sigma_ilr = sigma, order = nodes)[[2]]

X = seq(-0.7, 1, 0.01)
Y = seq(-0.5, 0.5, 0.01)
z2 = z1 = matrix(0, nrow = length(X), ncol = length(Y))
P0 = dnm(x, mu, sigma, nodes)
p = inv_ilr_coordinates(t(mu_max))
sigma_max = diag(1,3) / (mean(diag(sum(x) * (diag(c(p)) - t(p) %*% p))))
for(i in seq_along(X)){
  for(j in seq_along(Y)){
    z1[i,j] = exp(mvtnorm::dmvnorm(cbind(X[i],Y[j],-0.016761992), mu, sigma, log = T) +
                   lpmultinomial(x,
                                 inv_ilr_coordinates(cbind(X[i],Y[j],-0.016761992)),
                                 lconst = lpmultinomial_const(x))) / P0
    z2[i,j] = mvtnorm::dmvnorm(cbind(X[i],Y[j],-0.016761992), mu_max, sigma_max )
  }
}

XRANGE = range(pretty(X))
YRANGE = range(pretty(Y))
contour(X, Y, z1, xlim = XRANGE, ylim = YRANGE)
points(mu_exp[1], mu_exp[2], col='red')
points(mu_max[1], mu_max[2], col='blue')
par(new=TRUE)
contour(X, Y, z2, xlim = XRANGE, ylim = YRANGE, col = 'blue')
