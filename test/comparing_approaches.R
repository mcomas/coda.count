library(coda.base)

set.seed(1)
D = 3
N = 100
X = rlrnm(n = N, size = 20, c(0,0), diag(D-1))
B = ilr_basis(D)
#fit_0 = lapply(c(10:15), function(ord) c_fit_lrnm_hermite(X, B, order = ord))
fit_1 = lapply(c(10:15), function(ord) c_fit_lrnm_hermite_precision(X, B, order = ord, eps = 1e-08,
                                                                    em_max_steps = 30))
#fit_2 = lapply(c(10:15,100), function(ord) c_fit_lm_lrnm_hermite(X, B, matrix(1, ncol = 1, nrow = N), order = ord))
fit_3 = lapply(c(10:15), function(ord) c_fit_lm_lrnm_hermite_centered(X, B, matrix(1, ncol = 1, nrow = N),
                                                                      order = ord, eps = 1e-07,
                                                                      em_max_steps = 30))
fit_3_b = lapply(c(10:15), function(ord) c_fit_lm_lrnm_hermite_centered(X, B, matrix(1, ncol = 1, nrow = N),
                                                                        order = ord, eps = 1e-08,
                                                                        em_max_steps = 30))


#mu0 = t(sapply(fit_0, function(f) f[[1]]))
mu1 = t(sapply(fit_1, function(f) f[[1]]))
#mu2 = t(sapply(fit_2, function(f) f[[1]]))
mu3 = t(sapply(fit_3, function(f) f[[1]]))
mu3_b = t(sapply(fit_3_b, function(f) f[[1]]))
plot(mu1, xlim = range(c(mu1[,1], mu3[,1], mu3_b[,1])),
     ylim = range(c(mu1[,2], mu3[,2], mu3_b[,2])))
#points(mu2, col = 'blue')
#points(mu2, col = 'red')
points(mu3, col = 'blue')
points(mu3_b, col = 'yellow')
