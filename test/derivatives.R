library(coda.base)
x = c(0,3, 100)
mu = c(-200,0)
sigma = matrix(c(1,0.8,0.8,1), nrow = 2)
inv_sigma = solve(sigma)
B = alr_basis(3)
Binv = t(MASS::ginv(B))

I = 2
h = c(2,0.5)
deriv_1 = function(I, h, x, mu, inv_sigma, Binv){
  eBh = exp(Binv %*% h)
  w = sum(eBh)
  wBi = sum(Binv[,I] * eBh)
  -sum((h-mu) * inv_sigma[,I]) + sum(x * Binv[,I]) - sum(x) * wBi / w
}
deriv_1(1, h, x, mu, inv_sigma, Binv)
coda.count::lpnm_join_deriv(0, h, mu, inv_sigma, x)

coda.count::lpnm_join_deriv(1, h, mu, inv_sigma, x)
l_lrnm_join_d1(h, x, mu, inv_sigma, Binv)




J = 2
deriv_2 = function(I, J, h, x, mu, inv_sigma, Binv){
  eBh = exp(Binv %*% h)
  w = sum(eBh)
  #wBm = t(Binv) %*% eBh
  wBi = sum(Binv[,I] * eBh) # wBm[I]
  wBj = sum(Binv[,J] * eBh) # wBm[J]
  wBij = sum(Binv[,I] * Binv[,J] * eBh)
  -inv_sigma[I,J] - sum(x) * ( -wBi * wBj/ w^2 + wBij / w)
}
deriv_2(1, 1, h, x, mu, inv_sigma, Binv)
coda.count::lpnm_join_deriv2(0, 0, h, mu, inv_sigma, x)

l_lrnm_join_d2(h, x, mu, inv_sigma, Binv)
microbenchmark::microbenchmark(
  coda.count::lpnm_join_deriv2(0, 0, h, mu, inv_sigma, x),
  l_lrnm_join_d2(h, x, mu, inv_sigma, Binv)
)

