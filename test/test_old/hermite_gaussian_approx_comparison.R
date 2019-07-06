library(coda.count)
library(coda.base)

MU = c(1)
SIGMA = matrix(c(2), ncol = 1)
x = c(30900,0)
B = alr_basis(2)
Binv = t(MASS::ginv(B))

M = coda.count::c_obtain_moments_lrnm_hermite(matrix(x, nrow = 1), MU, SIGMA, B, 100)

# m1_hermite = c_m1_lrnm_hermite(x, MU, SIGMA, Binv, 100) / c_d_lrnm_hermite(x, MU, SIGMA, Binv, 100)
# m1_approx = c_posterior_alr_approximation(x, MU, SIGMA)[2]
#
# m2_hermite = c_m2_lrnm_hermite(x, MU, SIGMA, Binv, 100) / c_d_lrnm_hermite(x, MU, SIGMA, Binv, 100)
# m2c_hermite = m2_hermite - m1_hermite^2
# m2c_approx = c_posterior_alr_approximation(x, MU, SIGMA)[1]

# m1_hermite
# m1_approx
#
# m2c_hermite
# m2c_approx

to_ilr = function(h) coordinates(composition(h, 'alr'))



source('test/approximations.R')
#source('test/multinomial_to_mvnormal.R')
plot_approx(x, MU, SIGMA)
abline(v = to_ilr(m1_hermite), col = 'blue')
abline(v = to_ilr(m1_approx), col = 'red')

