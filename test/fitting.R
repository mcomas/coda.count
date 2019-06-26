library(coda.count)
library(coda.base)

X = rlrnm(n = 100, 20, c(0,0), diag(2))
B = alr_basis(3)

c_dm_fit_alpha(X)
c_fit_lrnm_hermite(X, B, order = 50)
