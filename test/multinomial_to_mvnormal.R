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
  f1 = sapply(h, function(h) gtools::ddirichlet(composition(h, 'alr'), ALPHA) * prod(composition(h, 'alr')))
  f2 = sapply(h, function(h) dnorm(h, N$MU, sqrt(N$SIGMA)))
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
  ALPHA = c(10,35)
  N = dirichlet_gaussian_approx_ILR(ALPHA)
  h = seq(-10, 5, 0.02)
  f1 = sapply(h, function(h) gtools::ddirichlet(composition(h), ALPHA) * sqrt(2) * prod(composition(h)))
  f2 = sapply(h, function(h) dnorm(h, N$MU, sqrt(N$SIGMA)))
  plot(h,f1, type = 'l', ylim = range(c(f1,f2)),xlim=c(-3,0), col='blue')
  points(h,f2, type='l', col='red')
}
gaussian_product = function(MU1, SIGMA1, MU2, SIGMA2){
  invSIGMA12 = solve(SIGMA1+SIGMA2)
  SIGMA3 = SIGMA1 %*% invSIGMA12 %*% SIGMA2
  MU3 = SIGMA2 %*% invSIGMA12 %*% t(t(MU1)) + SIGMA1 %*% invSIGMA12 %*% t(t(MU2))
  list(MU = MU3, SIGMA = SIGMA3)
}

plot_approx = function(X, MU, SIGMA){
  D = length(X)
  N1 = dirichlet_gaussian_approx(ALPHA = X + 1)
  MU_ALR = coordinates(composition(MU), 'alr')

  ISIGMA = solve(SIGMA)
  ILR2ALR = t(ilr_basis(length(X))) %*% alr_basis(length(X))
  SIGMA_ALR = t(ILR2ALR) %*% SIGMA %*% ILR2ALR

  APPROX = gaussian_product(N1$MU, N1$SIGMA, MU, SIGMA)
  prob = c_dlrnm_hermite_(X, MU, SIGMA, order = 1000)

  M1 = c_m1_hermite(x = X, mu = MU, sigma = SIGMA, order = 100)/prob
  M2 = c_m2_hermite(x = X, mu = MU, sigma = SIGMA, order = 100)/prob
  #list('MU' = M1, 'SIGMA' = M2 - M1^2)

  mx_alr = lpnm_join_maximum_alr(X, MU_ALR, solve(SIGMA_ALR), a = 0, eps = 0.00001, max_iter = 1000)
  mx = coordinates(composition(mx_alr, 'alr'))


  h = seq(-10, 10, 0.02)
  f = sapply(h, function(h) exp(lpnm_join(X, MU, ISIGMA, matrix(composition(h), nrow=1), matrix(h,nrow=1)))/prob)
  f1 = sapply(h, function(h) dnorm(h, MU, sqrt(SIGMA)))

  f2 = sapply(h, function(h) dnorm(h, APPROX$MU, sqrt(APPROX$SIGMA)))
  f3 = sapply(h, function(h) gtools::ddirichlet(composition(h), X+1) * sqrt(D) * prod(composition(h)))
  f4 = sapply(h, function(h) dnorm(h, N1$MU, sqrt(N1$SIGMA)))
  y_lim = range(c(f,f1,f2,f3,f4))
  f_max = pmax(f,f1,f2,f3,f4)
  x_lim = range(h[f_max>0.01])
  plot(h,f, type = 'l', ylim = y_lim, xlim = x_lim)
  abline(v = M1, col = 'black')
  #abline(v = APPROX$MU, col = 'green')
  #abline(v = mx, col = 'red')
  points(h, f1, type = 'l', col = 'blue')
  points(h, f2, type = 'l', col = 'green')
  points(h, f3, type = 'l', col = 'brown', lty = 2)
  points(h, f4, type = 'l', col = 'purple', lty = 2)
  legend('topright',
         legend = c('Prior', 'Posterior', 'Posterior approx.', 'Dirichlet', 'Dirichlet approx.'),
         col = c('blue', 'black', 'green', 'brown', 'purple'), cex=0.75, bty='n', lty=c(1,1,1,2,2))
}
plot_approx(X = c(0, 1), MU = c(2), SIGMA = diag(10, ncol = 1))

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
