library(coda.count)
library(coda.base)
load('testing/parliament2015.RData')
X = as.matrix(parliament2015[,3:7])
head(X)

mu0 = c(0,0)
sigma0 = matrix(c(1, -0.8, -0.8, 2), nrow = 2)

X = rlrnm(100, 50, mu0, sigma0)

Z = randtoolbox::sobol(n = 1000, dim = 2, normal = TRUE)
B = ilr_basis(ncol(X))

fit1 = fit_lrnm(X, B = B, method = 'montecarlo', Z = Z)
fit2 = fit_lrnm(X, B = B, method = 'hermite')
fit3 = fit_lrnm(X, B = B, method = 'laplace', max_iter = 10000)
fit4 = fit_mv_conditional_lrnm(X, B = B)
fit5 = fit_conditional_lrnm(X, B = B)
fit6 = fit_lrnm_vem_original(X)


Hfull = coordinates(X[apply(coordinates(X, B),1, function(x) all(is.finite(x))),])
colMeans(Hfull)

fit1[[1]]
fit2[[1]]
fit3[[1]]
fit4[[1]] %*% t(fit4$B1sub)
fit5[[1]] %*% t(fit5$B1sub)
coordinates(composition(fit6[[1]], 'alr'), B)

colMeans(coordinates(fit6$P))

sum(diag(cov(Hfull)))
sum(diag(fit1[[2]]))
sum(eigen(fit2[[2]])$values)
sum(eigen(fit3[[2]])$values)
eigen(fit4$B1sub %*% fit4[[2]] %*% t(fit4$B1sub))$values
eigen(fit5$B1sub %*% fit5[[2]] %*% t(fit5$B1sub))$values
eigen(cov(coordinates(fit6$P, B)))$values
