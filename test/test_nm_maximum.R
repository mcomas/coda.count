library(Rcpp)
sourceCpp('src/nm_maximum.cpp')
library(microbenchmark)

k = 50
x = rep(1, k+1)
mu = rep(0, k)
sigma = diag(k)
microbenchmark(
  example_NewtonDescentSolver(x, mu, sigma),
  bfgsSolver(x, mu, sigma))
