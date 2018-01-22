set.seed(1)
library(coda.count)

MU.orig = c(0,0)
SIGMA.orig = diag(2)
X = rlrnm(MU.orig, SIGMA.orig, 100, 10)

library(Rcpp)
sourceCpp('src/nm_maximum.cpp')

vv = c_lrnm_fit_maximum_alr(X, MU.orig, SIGMA.orig)
