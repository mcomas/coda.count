library(coda.count)
library(coda.base)
library(dplyr)
library(tidyr)
library(ggtern)


lrnm_posterior = function(MU, SIGMA, x){
  X_MU = composition(MU)
  ILR2ALR = t(ilr_basis(3)) %*% alr_basis(3)
  ALR2ILR = solve(ILR2ALR)

  MU_ALR = MU %*% ILR2ALR
  #composition(MU_ALR, 'alr')
  #composition(MU, 'ilr')

  SIGMA_ALR = t(ILR2ALR) %*% SIGMA %*% ILR2ALR
  a_opt = coda.count::lpnm_join_maximum_alr(x, MU_ALR, solve(SIGMA_ALR),
                                            a = MU_ALR, eps = 0.0001, max_iter = 100)
  x_opt = composition(a_opt[,1], 'alr')
  #ldnormal(matrix(h, nrow = 1), MU, solve(SIGMA))
  #mvtnorm::dmvnorm(matrix(h, nrow = 1), mean = MU, sigma = SIGMA, log = TRUE)

  Pc = simplex_lattice(100, 3)
  Pc = Pc[rowSums(Pc == 0) == 0,]
  P = Pc/rowSums(Pc)
  H = coordinates(P)


  df = bind_cols(as_tibble(P, .name_repair = ~paste0('p',1:3)),
                 as_tibble(H, .name_repair = ~paste0('h',1:2)))
  df$f_prior = scale(exp(ldnormal(H, MU, solve(SIGMA))))[,1]
  df$f_posterior = scale(exp(lpnm_join(x, MU, solve(SIGMA), P, H)) / dlrnm(x, MU, SIGMA))[,1]

  dfplot =  df %>%
    gather(distribution, density, starts_with('f_')) %>%
    mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior'),
                                 labels = c('Prior', 'Posterior')))
  dopt = bind_rows(
    tibble(x = x_opt[2], y = x_opt[1], z = x_opt[3], distribution = 'f_posterior'),
    tibble(x = X_MU[2], y = X_MU[1], z = X_MU[3], distribution = 'f_prior')) %>%
    mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior'),
                                 labels = c('Prior', 'Posterior')))
  list(data = dfplot, opt = dopt)
}

MU = c(0,0)
SIGMA = matrix(c(1,-0.8,
                 -0.8,1), nrow = 2) * 0.25
x = c(10,1,20)

res = lrnm_posterior(MU, SIGMA, x)

ggtern() +
  geom_point(data=res$data, aes(x=p2,y=p1,z=p3, col=density), size = 2) +
  geom_mask() +
  geom_point(data = data.frame(x = x[2], y = x[1], z = x[3]), aes(x=x,y=y,z=z), col = 'black') +
  geom_point(data = res$opt,
             aes(x=x,y=y,z=z), col = 'blue') +
  scale_color_continuous(low = 'white', high = 'red') +
  facet_wrap(~distribution) + theme(legend.position = 'none')

if(FALSE){
  h = coordinates(x)
  ggplot() +
    geom_point(data=dfplot, aes(x=h1,y=h2, col=density), size = 1) +
    geom_point(data = data.frame(x = h[1], y = h[2]), aes(x=x,y=y), col = 'black') +
    scale_color_continuous(low = 'white', high = 'red') +
    facet_wrap(~distribution) +
    coord_equal() +
    theme_minimal() + theme(legend.position = 'none')
}
###########

# X = apply(simplex_lattice(10,2), 2, rev)
# sum(p <- apply(X, 1, dlrnm, -0.5, 0.1))
# names(p) = sprintf("(%d,%d)", X[,1], X[,2])
# barplot(p, cex.axis = 0.8, cex.names = 0.8)
