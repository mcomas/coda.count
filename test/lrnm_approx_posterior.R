posterior_approx = function(ALPHA, x){
  ALR2ILR = MASS::ginv(alr_basis(3)) %*% ilr_basis(3)

  MU_ALR = digamma(ALPHA[1:2])-digamma(ALPHA[3])
  MU = coordinates(composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR

  SIGMA_ALR = diag(trigamma(ALPHA[1:2])) + trigamma(ALPHA[3])
  SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR

  h = coordinates(x)

  ldnormal(matrix(h, nrow = 1), MU, solve(SIGMA))
  mvtnorm::dmvnorm(matrix(h, nrow = 1), mean = MU, sigma = SIGMA, log = TRUE)

  Pc = simplex_lattice(100, 3)
  Pc = Pc[rowSums(Pc == 0) == 0,]
  P = Pc/rowSums(Pc)
  H = coordinates(P)


  df = bind_cols(as_tibble(P, .name_repair = ~paste0('p',1:3)),
                 as_tibble(H, .name_repair = ~paste0('h',1:2)))
  df1 = df
  df1$f_prior = scale(exp(ldnormal(H, MU, solve(SIGMA))))[,1]
  df1$f_posterior = scale(exp(lpnm_join(x, MU, solve(SIGMA), P, H)) / dlrnm(x, MU, SIGMA))[,1]

  df2 = df
  df2$f_prior = scale(apply(P, 1, ddirichlet, ALPHA))[,1]
  df2$f_posterior = scale(apply(P, 1, ddirichlet, ALPHA + x))[,1]

  dfplot =  bind_rows(df1 %>% mutate(type = 'Dirichlet'),
                      df2 %>% mutate(type = 'Log-Ratio-Normal')) %>%
    gather(distribution, density, starts_with('f_')) %>%
    mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior'),
                                 labels = c('Prior', 'Posterior')))
  list(data = dfplot, m = composition(MU))
}

X = c(0,20,40)
dfplot = posterior_approx(ALPHA = 5*c(1,3,2), x = X)
ggtern() +
  geom_point(data=dfplot$data, aes(x=p2,y=p1,z=p3, col=density), size = 2) +
  geom_mask() +
  geom_point(data = data.frame(x = X[2], y = X[1], z = X[3]),
             aes(x=x,y=y,z=z), col = 'black') +
  geom_point(data = data.frame(x = dfplot$m[2], y = dfplot$m[1], z = dfplot$m[3]),
             aes(x=x,y=y,z=z), col = 'blue') +
  scale_color_continuous(low = 'white', high = 'red') +
  facet_grid(type~distribution) + theme_bw() +
  theme(legend.position = 'none')

if(FALSE){
  H = coordinates(X)
  ggplot() +
    geom_point(data=dfplot, aes(x=h1,y=h2, col=density), size = 1) +
    geom_point(data = data.frame(x = H[1], y = H[2]), aes(x=x,y=y), col = 'black') +
    scale_color_continuous(low = 'white', high = 'red') +
    facet_grid(type~distribution) +
    coord_equal() +
    theme_void() + theme(legend.position = 'none')
}
