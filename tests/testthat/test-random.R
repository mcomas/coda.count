
test_that("multinomial", {
  X1 = rmultinomial(10, 100, c(0.2,0.3,0.5))
  expect_that(X1, is_a("matrix"))
  expect_type(X1, 'integer')
  probs = matrix(c(0.2, 0.4, 0.4,
                   0.5, 0.1, 0.4), byrow = TRUE, nrow = 6, ncol = 3)
  size = c(10, 500, 10000)
  X2 = rmultinomial(size = size, p = probs)
  expect_that(X2, is_a("matrix"))
  expect_type(X2, 'integer')
})

test_that("dirichlet-multinomial", {
  X1 = rdm(n = 100, 10, c(1.5, 2, 2.5))
  expect_that(X1, is_a("matrix"))
  #expect_type(X1, 'integer')
  probs = matrix(c(0.2, 0.4, 0.4,
                   0.5, 0.1, 0.4), byrow = TRUE, nrow = 6, ncol = 3)
  size = c(10, 500, 10000)
  X2 = rdm(size = c(1000, 100, 10, 2, 1), alpha = c(1,1,1))
  expect_that(X2, is_a("matrix"))
  #expect_type(X2, 'integer')
})
