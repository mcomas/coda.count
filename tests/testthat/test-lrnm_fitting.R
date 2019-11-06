test_that("fitting", {
  set.seed(1)
  X = rlrnm(n = 10, size = 100000, mu = c(0,0), sigma = diag(2))
  expect_type(fit_lrnm(X, method = 'montecarlo'), "list")
})
