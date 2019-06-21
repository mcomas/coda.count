library(coda.base)
library(coda.count)
dirichlet_gaussian_approx = function(ALPHA){
  D = length(ALPHA)
  ILR2ALR = t(ilr_basis(D)) %*% alr_basis(D)
  ALR2ILR = solve(ILR2ALR)

  MU_ALR = digamma(ALPHA[-D]) - digamma(ALPHA[D])
  MU = coordinates(composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR

  SIGMA_ALR = diag(trigamma(ALPHA[-D]), nrow = D-1) + trigamma(ALPHA[D])
  SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR
  return(list(MU = MU, SIGMA = SIGMA))
}
gaussian_product = function(MU1, SIGMA1, MU2, SIGMA2){
  invSIGMA12 = solve(SIGMA1+SIGMA2)
  SIGMA3 = SIGMA1 %*% invSIGMA12 %*% SIGMA2
  MU3 = SIGMA2 %*% invSIGMA12 %*% t(t(MU1)) + SIGMA1 %*% invSIGMA12 %*% t(t(MU2))
  list(MU = MU3, SIGMA = SIGMA3)
}
gaussian_division = function(MU2, SIGMA2, MU3, SIGMA3){
  invSIGMA3 = solve(SIGMA3)
  invSIGMA2 = solve(SIGMA2)
  SIGMA1 = solve(invSIGMA3 - invSIGMA2)
  MU1 = SIGMA1 %*% invSIGMA3 %*% MU3 - SIGMA1 %*% invSIGMA2 %*% MU2
  list(MU = MU1, SIGMA = SIGMA1)
}

if(FALSE){
  dirichlet_gaussian_approx_ALR = function(ALPHA){
    D = length(ALPHA)

    MU_ALR = digamma(ALPHA[-D]) - digamma(ALPHA[D])
    SIGMA_ALR = diag(trigamma(ALPHA[-D]), nrow = D-1) + trigamma(ALPHA[D])

    return(list(MU = MU_ALR, SIGMA = SIGMA_ALR))
  }
  ALPHA = c(10,35)
  N = dirichlet_gaussian_approx_ALR(ALPHA)
  h = seq(-10, 5, 0.02)
  f1 = sapply(h, function(h) gtools::ddirichlet(composition(h, 'alr'), ALPHA))
  f2 = sapply(h, function(h) dnorm(h, N$MU, sqrt(N$SIGMA)) / prod(composition(h, 'alr')))
  plot(h,f1, type = 'l', ylim = range(c(f1,f2)),xlim=c(-3,0), col='blue')
  points(h,f2, type='l', col='red')
}
if(FALSE){
  dirichlet_gaussian_approx_ILR = function(ALPHA){
    D = length(ALPHA)
    ILR2ALR = t(ilr_basis(D)) %*% alr_basis(D)
    ALR2ILR = solve(ILR2ALR)

    MU_ALR = digamma(ALPHA[-D]) - digamma(ALPHA[D])
    MU = coordinates(composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR

    SIGMA_ALR = diag(trigamma(ALPHA[-D]), nrow = D-1) + trigamma(ALPHA[D])
    SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR
    return(list(MU = MU, SIGMA = SIGMA))
  }
  ALPHA = c(10,2)
  N = dirichlet_gaussian_approx_ILR(ALPHA)
  h = seq(-2, 5, 0.02)
  f1_ = function(h) sapply(h, function(h) gtools::ddirichlet(composition(h), ALPHA) )
  f2_ = function(h) sapply(h, function(h) dnorm(h, N$MU, sqrt(N$SIGMA)) / (sqrt(2) * prod(composition(h))))
  f1 = f1_(h) / integrate(f1_, lower = -5, upper = 5)$value
  f2 = f2_(h) / integrate(f2_, lower = -5, upper = 5)$value
  MU2 = 0
  SIGMA2 = 1/(2 * pi * prod(composition(0))^2)
  N1 = gaussian_division(MU2, SIGMA2, N$MU, N$SIGMA)
  f3 = sapply(h, function(h) dnorm(h, N1$MU, sqrt(N1$SIGMA)))
  plot(h,f1, type = 'l', ylim = range(c(f1,f2)),xlim=c(-2,5), col='blue')
  points(h,f2, type='l', col='red')
  points(h,f3, type='l', col='green')
}

plot_approx = function(X, MU, SIGMA){
  if(FALSE){
    X = c(10, 2)
    MU = c(4)
    SIGMA = diag(1, ncol = 1)
  }
  D = length(X)
  dlrnm = c_dlrnm_hermite_(X, MU, SIGMA, order = 1000)
  ISIGMA = solve(SIGMA)
  f_posterior_ = function(h_) sapply(h_, function(h) exp(lpnm_join(X, MU, ISIGMA, matrix(composition(h), nrow=1), matrix(h,nrow=1)))/dlrnm)
  f_multinomial_ = function(h_, x) sapply(h_, function(h) gtools::ddirichlet(composition(h), x + 1) )
  h = seq(MU-30, MU+30, 0.02)

  N_dirichlet_lebesgue = dirichlet_gaussian_approx(X+1)

  MU_logistic = 0
  SIGMA_logistic = 1/(2 * pi * prod(composition(0))^2)
  N_dirichlet_aitchison = gaussian_division(MU_logistic, SIGMA_logistic,
                                            N_dirichlet_lebesgue$MU, N_dirichlet_lebesgue$SIGMA)
  N_posterior_approx = gaussian_product(MU, SIGMA,
                                        N_dirichlet_aitchison$MU, N_dirichlet_aitchison$SIGMA)


  f_posterior = f_posterior_(h)
  f_prior = dnorm(h, MU, sqrt(SIGMA))
  f_multinomial = f_multinomial_(h, X) / integrate(f_multinomial_, -10, 10, x = X)$value
  f_dirichlet_lebesgue = dnorm(h, N_dirichlet_lebesgue$MU, sqrt(N_dirichlet_lebesgue$SIGMA))
  f_dirichlet_aitchison = dnorm(h, N_dirichlet_aitchison$MU, sqrt(N_dirichlet_aitchison$SIGMA))
  f_posterior_approximation = dnorm(h, N_posterior_approx$MU, sqrt(N_posterior_approx$SIGMA))

  f_max = pmax(f_posterior, f_prior, f_multinomial, f_dirichlet_aitchison, f_posterior_approximation)
  x_lim = range(h[f_max>0.01])
  y_lim = range(c(f_posterior, f_prior, f_multinomial, f_dirichlet_aitchison, f_posterior_approximation))
  plot(h,f_posterior, type = 'l', ylim = y_lim, xlim = x_lim, lwd = 3)
  points(h, f_prior, type = 'l', col = 'blue', lty = 2)
  points(h, f_multinomial, type = 'l', col = 'red')
  points(h, f_dirichlet_aitchison, type = 'l', col = 'red', lty = 2)
  #points(h, f_dirichlet_lebesgue, type = 'l', col = 'red', lty = 3)
  points(h, f_posterior_approximation, type = 'l', col = 'green', lty = 2)
  legend('topright',
         legend = c('Posterior', 'Posterior approx.', 'Prior',
                    'Multinomial', 'Multinomial approx.'),
         col = c('black', 'green', 'blue', 'red', 'red'), cex=0.75, bty='n', lty=c(1,2,1,1,2))

}

plot_approx(X = c(1, 2000), MU = c(2), SIGMA = diag(1, ncol = 1))
plot_approx(X = c(200, 1), MU = c(5), SIGMA = diag(1, ncol = 1))
plot_approx(X = c(20, 0), MU = c(2), SIGMA = diag(1, ncol = 1))
plot_approx(X = c(10,10), MU = c(2), SIGMA = diag(1, ncol = 1))

## Bad parameters
# If Dirichlet is almost constant within prior domain, posterior remains as the prior.
# X = 1000*c(100, 0); MU = c(20); SIGMA = diag(1, ncol = 1)
#
# X = c(100, 2); MU = c(10); SIGMA = diag(1, ncol = 1)
# plot_approx(X = c(100, 1), MU = c(10), SIGMA = diag(1, ncol = 1))
if(FALSE){
  library(manipulate)
  manipulate(
    plot_approx(X = a * c(10, 0), MU = c(2), SIGMA = diag(1, ncol = 1)),
    a = slider(1, 1000, step = 100)
  )
}
if(FALSE){
  h = seq(-5, 5, 0.02)
  logistic_ilr = function(h) sapply(h, function(h_) prod(composition(h_)) * sqrt(2))
  logistic_alr = function(h) sapply(h, function(h_) prod(composition(h_, 'alr')))
  integrate(logistic_ilr, -100, 100)$value
  integrate(logistic_alr, -100, 100)$value
  f_ilr = logistic_ilr(h)
  f_ilr_approx = dnorm(h, mean = 0, sd = 1/sqrt(2 * pi * logistic_ilr(0)^2))
  f_alr = logistic_alr(h)
  f_alr_approx = dnorm(h, mean = 0, sd = 1.6)
  plot(h,f_ilr, type = 'l')
  points(h, f_ilr_approx, type = 'l', col = 'green')
  points(h, f_alr, type = 'l', col = 'red')
  points(h, f_alr_approx, type = 'l', col = 'blue')
}
if(FALSE){
  library(ggplot2)
  H = expand.grid(h1 = seq(-2, 2, 0.05), h2 = seq(-2, 2, 0.05))
  logistic_ilr = function(H) apply(H, 1, function(h_) prod(composition(h_)) * sqrt(3))
  H$f = logistic_ilr(H)
  hmax = logistic_ilr(matrix(0,ncol=2))
  H$f_gaussian = mvtnorm::dmvnorm(H[,1:2], mean = c(0,0), diag(1/(hmax*2*pi), 2))
  ggplot(data=H) +
    geom_point(aes(x=h1, y = h2, col = f))
  ggplot(data=H) +
    geom_point(aes(x=h1, y = h2, col = f_gaussian))
}
#
# p = 0.001
# n = 10
# x = 0:n
# y = n-x
# xc = seq(0,n,0.001)
# yc = n-xc
#
# fx = dbinom(x, n, p)
# fx_n = dnorm(xc, n * p, sqrt(n * p * (1-p)))
# plot(x,fx, ylim = range(c(fx,fx_n)))
# points(xc,fx_n, col = 'red', type='l')
#
# fy = dbinom(y, n, 1-p)
# fy_n = dnorm(yc, n * (1-p), sqrt(n * p * (1-p)))
# plot(y,fy, ylim = range(c(fy,fy_n)))
# points(yc,fy_n, col = 'red', type='l')
