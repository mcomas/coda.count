library(coda.count)
library(coda.base)
library(microbenchmark)
D = 2
B = ilr_basis(D)

X = rlrnm(n = 100, 10000, rep(0,D-1), diag(D-1), B = B)

c_dm_fit_alpha(X)
c_fit_lrnm_hermite(X, B, order = 30, eps = 0.000001)
c_fit_lrnm_hermite_precision(X, B, order = 15, eps = 0.000001)

# microbenchmark(
#   c_fit_lrnm_hermite(X, B, order = 5),
#   c_fit_lrnm_hermite_fast(X, B, order = 5))
# Unit: milliseconds
#                                     expr      min       lq     mean   median       uq      max neval
#      c_fit_lrnm_hermite(X, B, order = 5) 22.49687 22.85240 23.28217 23.00135 23.19763 28.67885   100
# c_fit_lrnm_hermite_fast(X, B, order = 5) 10.59516 10.68236 10.90090 10.74042 10.85607 14.14484   100
c_fit_lrnm_gaussian_approx(X, B)
