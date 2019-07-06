
plot_ = function(mu1, mu2, var1, var2, corr, I){
  SIZE = 100
  mu = c(mu1,mu2)
  sigma = matrix(c(var1,corr*sqrt(var1)*sqrt(var2),
                   corr*sqrt(var1)*sqrt(var2),var2), nrow = 2)
  SL3 = simplex_lattice(SIZE, 3)
  lM = lapply(split(SL3, 1:nrow(SL3)), c_posterior_approximation, mu, sigma, t(MASS::ginv(B)))
  H_MAX = as.data.frame(t(sapply(lM, function(l)l[,3])))
  MAX = composition(H_MAX, B)
  H_E = as.data.frame(t(c_obtain_moments_lrnm_hermite(SL3, mu, sigma, ilr_basis(3), order = 10)[[1]]))
  E = composition(H_E, B)
  p = ggplot() +
    geom_point(data=H_MAX,  aes(x = V1, y = V2), col = 'grey', size = 0.2)+
    geom_point(data=H_E,  aes(x = V1, y = V2), size = 0.2)
  for(i in I){
    p = p +
      geom_point(data=H_E[i,,drop=FALSE],  aes(x = V1, y = V2)) +
      geom_point(data=H_MAX[i,,drop=FALSE],  aes(x = V1, y = V2), col = 'red') +
      geom_path(data = ellipse(lM[[i]][,3], lM[[i]][1:2,1:2], 0.95),
                aes(x = V1, y = V2), alpha = 0.5, col = 'red')
  }
  p +
    geom_point(aes(x = mu[1], y = mu[2]), col = 'blue') +
    geom_path(data = ellipse(mu, sigma, 0.95), aes(x = V1, y = V2), col = 'blue') +
    coord_equal()
}

plot_(0,-2,1,1,0, 1) +
  coord_cartesian(xlim = c(-1,5), ylim = c(-1,1))





# ```
#
# ****
#
#   ```{r, fig.width=6, fig.height=4.8, fig.align='center'}
# SL3 = simplex_lattice(SIZE, 3)
# M = c_obtain_moments_lrnm_hermite(SL3, mu, sigma, ilr_basis(3), order = 10)
#
# Binv = t(MASS::ginv(B))
# lM = lapply(split(SL3, 1:nrow(SL3)), c_posterior_approximation, mu, sigma, Binv)
#
# # H1 = as.data.frame(t(M[[1]]))
# # M1 = composition(H1)
# ellipse = function(mu, sigma, p){
#   s = -2 * log(1 - p);
#   ev = eigen(sigma * s)
#   t_ = seq(0, 2 * pi, length.out = 500)
#   a = mu + t((t(ev$vectors) * sqrt(ev$values))) %*% rbind(cos(t_), sin(t_))
#   as.data.frame(t(a))
# }
# coda_ellipse = function(mu, sigma, p, B = ilr_basis(length(mu) + 1)){
#   composition(ellipse(mu, sigma, p), B)
# }
#
# ggtern() +
#   geom_point(data = as.data.frame(rlrnormal(100, mu, sigma)),
#              aes(x = V1, y = V2, z = V3), alpha = 0.5) +
#   geom_path(data = coda_ellipse(mu, sigma, 0.95), aes(x = x1, y = x2, z = x3))
# ```
#
# ```{r}
# p = ggtern() +
#   geom_point(aes(x = 0.4359516, y = 0.4359516, z = 0.1280968), col = 'blue') +
#   geom_path(data = coda_ellipse(mu, sigma, 0.95), aes(x = x1, y = x2, z = x3), col = 'blue')
# for(i in 1){
#   p = p +
#     geom_point(data = M1[i,], aes(x = x1, y = x2, z = x3), alpha = 0.5) +
#     geom_path(data = coda_ellipse(lM[[i]][,3], lM[[i]][1:2,1:2], 0.95), aes(x = x1, y = x2, z = x3), alpha = 0.5)
# }
# p
# ```
#
# ```{r}
# p = ggplot() +
#   geom_point(aes(x = mu[1], y = mu[2]), col = 'blue') +
#   geom_path(data = ellipse(mu, sigma, 0.95), aes(x = V1, y = V2), col = 'blue')
# x = c_posterior_approximation(c(10,2,1), mu, sigma, Binv)
# p = p +
#   geom_point(data = H1[i,], aes(x = V1, y = V2), alpha = 0.5) +
#   geom_path(data = ellipse(lM[[i]][,3], lM[[i]][1:2,1:2], 0.95), aes(x = V1, y = V2), alpha = 0.5)
#
# p
#
# p = ggplot() +
#   geom_point(aes(x = mu[1], y = mu[2]), col = 'blue') +
#   geom_path(data = ellipse(mu, sigma, 0.95), aes(x = V1, y = V2), col = 'blue')
# for(i in 1){
#   p = p +
#     #geom_point(data = coordinates(SL3[i,]), aes(x = V1, y = V2), alpha = 0.5) +
#     geom_point(data = H1[i,], aes(x = V1, y = V2), alpha = 0.5) +
#     geom_path(data = ellipse(lM[[i]][,3], lM[[i]][1:2,1:2], 0.95), aes(x = V1, y = V2), alpha = 0.5)
# }
# p
# ```
