test_that("fitting", {
  set.seed(1)
  X = rlrnm(n = 5, size = 1000, mu = c(0,0), sigma = diag(2))
  fit_lrnm(X, method = 'montecarlo')
})
