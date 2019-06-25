x = c(0,3, 100)
mu = c(-2,0)
sigma = matrix(c(1,0.8,0.8,1), nrow = 2)
inv_sigma = solve(sigma)
B = alr_basis(3)
Binv = t(MASS::ginv(B))

I = 2
h = c(2,0)
deriv_1 = function(I, h, x, mu, inv_sigma, Binv){
  w = sum(exp(Binv %*% h))
  wBi = sum(Binv[,I] * exp(Binv %*% h))
  -sum((h-mu) * inv_sigma[,I]) + sum(x * Binv[,I]) - sum(x) * wBi / w
}
deriv_1(1, h, x, mu, inv_sigma, Binv)
coda.count::lpnm_join_deriv(0, h, mu, inv_sigma, x)

J = 2
deriv_2 = function(I, J, h, x, mu, inv_sigma, Binv){
  w = sum(exp(Binv %*% h))
  wBi = sum(Binv[,I] * exp(Binv %*% h))
  wBj = sum(Binv[,J] * exp(Binv %*% h))
  wBij = sum(Binv[,I] * Binv[,J] * exp(Binv %*% h))
  -inv_sigma[I,J] - sum(x) * ( -wBi * wBj/ w^2 + wBij / w)
}
deriv_2(1, 1, h, x, mu, inv_sigma, Binv)
coda.count::lpnm_join_deriv2(0, 0, h, mu, inv_sigma, x)
