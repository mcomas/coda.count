posterior_approx = function(ALPHA, x){

  ILR2ALR = t(ilr_basis(3)) %*% alr_basis(3)
  ALR2ILR = solve(ILR2ALR)

  MU_ALR = digamma(ALPHA[1:2])-digamma(ALPHA[3])
  MU = coordinates(X_MU <- composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR


  SIGMA_ALR = diag(trigamma(ALPHA[1:2])) + trigamma(ALPHA[3])
  SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR

  a_opt = coda.count::lpnm_join_maximum_alr(x, MU_ALR, solve(SIGMA_ALR),
                                            a = MU_ALR, eps = 0.000001, max_iter = 1000)
  x_opt = composition(a_opt[,1], 'alr')

  # ldnormal(matrix(h, nrow = 1), MU, solve(SIGMA))
  # mvtnorm::dmvnorm(matrix(h, nrow = 1), mean = MU, sigma = SIGMA, log = TRUE)

  Pc = simplex_lattice(120, 3)
  Pc = Pc[rowSums(Pc == 0) == 0,]
  P = Pc/rowSums(Pc)
  H = coordinates(P)

  df = bind_cols(as_tibble(P, .name_repair = ~paste0('p',1:3)),
                 as_tibble(H, .name_repair = ~paste0('h',1:2)))
  df1 = df
  df1$f_prior = scale(exp(ldnormal(H, MU, solve(SIGMA))))[,1]
  df1$f_posterior = scale(exp(lpnm_join(x, MU, solve(SIGMA), P, H)) / dlrnm(x, MU, SIGMA))[,1]

  df2 = df
  df2$f_prior = scale(apply(P, 1, gtools::ddirichlet, ALPHA))[,1]
  df2$f_posterior = scale(apply(P, 1, gtools::ddirichlet, ALPHA + x))[,1]

  dfplot =  bind_rows(df1 %>% mutate(type = 'Log-Ratio-Normal'),
                      df2 %>% mutate(type = 'Dirichlet')) %>%
    gather(distribution, density, starts_with('f_')) %>%
    mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior'),
                                 labels = c('Prior', 'Posterior')))

  xdirich = ALPHA
  xdirich2 = ALPHA + x
  dmax = bind_rows(
    tibble(x = X_MU[2], y = X_MU[1], z = X_MU[3], distribution = 'f_prior') %>%
      mutate(type = 'Log-Ratio-Normal'),
    tibble(x = xdirich[2], y = xdirich[1], z = xdirich[3], distribution = 'f_prior') %>%
      mutate(type = 'Dirichlet'),
    tibble(x = x_opt[2], y = x_opt[1], z = x_opt[3], distribution = 'f_posterior') %>%
      mutate(type = 'Log-Ratio-Normal'),
    tibble(x = xdirich2[2], y = xdirich2[1], z = xdirich2[3], distribution = 'f_posterior') %>%
      mutate(type = 'Dirichlet')) %>%
    mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior'),
                                 labels = c('Prior', 'Posterior')))
  dmax = bind_cols(dmax, dmax %>% select(y, x, z) %>% coordinates(label = 'h'))
  list(data = dfplot, m = dmax)
}

X = c(1,20,40)
dfplot = posterior_approx(ALPHA = 5*c(1,3,2), x = X)
ggtern() +
  geom_point(data=dfplot$data, aes(x=p2,y=p1,z=p3, col=density), size = 2) +
  geom_mask() +
  geom_point(data = data.frame(x = X[2], y = X[1], z = X[3]),
             aes(x=x,y=y,z=z), col = 'black') +
  geom_point(data = dfplot$m,
             aes(x=x,y=y,z=z), col = 'blue') +
  scale_color_continuous(low = 'white', high = 'red') +
  facet_grid(type~distribution) + theme_bw() +
  theme(legend.position = 'none')

if(FALSE){
  H = coordinates(X)
  ggplot() +
    geom_point(data=dfplot$data, aes(x=h1,y=h2, col=density), size = 1) +
    geom_point(data = dfplot$m,
               aes(x=h1,y=h2), col = 'blue') +
    geom_point(data = data.frame(x = H[1], y = H[2]), aes(x=x,y=y), col = 'black') +
    scale_color_continuous(low = 'white', high = 'red') +
    facet_grid(type~distribution) +
    coord_equal() +
    theme_void() + theme(legend.position = 'none')
}
