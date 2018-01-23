set.seed(1)
library(coda.count)
library(coda.base)

MU.orig = c(0,0)
SIGMA.orig = diag(2)
X = rlrnm(MU.orig, SIGMA.orig, 10, 100)

#library(Rcpp)
#sourceCpp('src/nm_maximum.cpp')

vv = c_lrnm_fit_maximum_alr(X, MU.orig, SIGMA.orig)
P = exp(cbind(vv,0))
P = P/rowSums(P)

colMeans(coordinates(P))
cov(coordinates(P))
