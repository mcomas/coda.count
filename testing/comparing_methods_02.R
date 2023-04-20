library(coda.count)
load('testing/parliament2015.RData')
X = as.matrix(parliament2015[,3:7])
head(X)

mu0 = c(0,0)
sigma0 = matrix(c(1, -0.8, -0.8, 2), nrow = 2)

compare_methods = function(mu,sigma,n,m){
  # mu = mu0
  # sigma = sigma0
  # n = 100
  # m = 10
  X = rlrnm(n, m, mu, sigma)

  fit1 = fit_lrnm(X, method = 'mc')
  fit2 = fit_lrnm(X, method = 'hermite')
  fit3 = fit_lrnm(X, method = 'laplace')
  fit4 = fit_conditional_lrnm(X, method = 'mc')
  fit5 = fit_conditional_lrnm(X, method = 'hermite')
  fit6 = fit_one_dimensional_conditional_lrnm(X)
  fit7 = fit_lrnm(X, method = 'vem')

  # fit6 = fit_vem_lrnm(X, B)
  H = coordinates(X)
  H[!is.finite(H)] = NA

  list('Full' = list('mu' = colMeans(H, na.rm = TRUE),
                     'sigma' = cov(H, use = 'complete.obs')),
       'lrnm-montecarlo' = fit1,
       'lrnm-hermite' = fit2,
       'lrnm-laplace' = fit3,
       'lrnm-cond-mc' = fit4,
       'lrnm-cond-hermite' = fit5,
       'lrnm-low-dim' = fit6,
       'lrnm-vem' = fit7)
}
lres = replicate(9, compare_methods(mu0, sigma0, 100, 10), simplify = FALSE)
lmu = lapply(seq_along(lres), function(i){
  dres = as.data.frame(t(sapply(lres[[i]], function(fit) fit$mu)))
  dres$method = names(lres[[i]])
  row.names(dres) = NULL
  dres
})
library(tidyverse)
dmu = bind_rows(lmu, .id = 'id')
ggplot(data=dmu) +
  geom_label(aes(x = ilr1, y = ilr2, label = id)) +
  facet_wrap(~method)
ggplot(data=dmu) +
  geom_point(aes(x = ilr1, y = ilr2, col= method)) +
  facet_wrap(~id, scales = 'free')

lsigma = lapply(seq_along(lres), function(i){
  dres = lapply(lres[[i]], function(fit){
    as.data.frame(ellipse::ellipse(fit$sigma))
  })
  names(dres) = names(lres[[i]])
  bind_rows(dres, .id = 'method')
})
dsigma = bind_rows(lsigma, .id = 'id')
ggplot(data=dsigma) +
  geom_path(aes(x = ilr1, y = ilr2, col = method)) +
  facet_wrap(~id)
