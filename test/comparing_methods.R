library(coda.count)
load('data/parliament2015.RData')
X = as.matrix(parliament2015[,3:7])
head(X)


X = rlrnm(100, 50, c(0,0), matrix(c(1, -0.8, -0.8, 2), nrow = 2))

Z = randtoolbox::sobol(n = 10000, dim = 2, normal = TRUE)
fit1 = fit_lrnm(X, method = 'montecarlo', Z = Z)
fit2 = fit_lrnm(X, method = 'hermite')
fit3 = fit_lrnm(X, method = 'laplace', max_iter = 10000)

fit1[1:2]
fit2[1:2]
fit3[1:2]
