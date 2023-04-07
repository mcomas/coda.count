F_optim = function(w, x, mu, sigma){
  I = 1:length(mu)
  m = w[+I]
  V = w[-I]
  invSigma = MASS::ginv(sigma)
  xi = sum(exp(m + 0.5 * V))
  -as.vector(sum(x[-length(x)] * m) - sum(x) * (xi^-1 * sum(exp(m + 0.5 * V)) - 1 + log(xi) ) +
              - 0.5 * (log(det(sigma)) + t(m-mu) %*% invSigma %*% (m-mu) + sum(diag(invSigma) * V)) +
              + 0.5 * (sum(log(V)) + length(m)))
}
gr_optim = function(w, x, mu, sigma){
  I = 1:length(mu)
  m = w[+I]
  V = w[-I]

  invSigma = MASS::ginv(sigma)

  gr1 = x[-length(x)] - as.vector(invSigma %*% (m-mu)) -
    as.vector(sum(x) * xi^-1 * exp(m + 0.5 * V))
  gr2 =   as.vector(1/V - V * diag(invSigma) - sum(x) * xi^-1 * exp(m + 0.5 * V^2) * V)
  -c(gr1,gr2)
}
#' @export
fit_vem_lrnm_optim = function(X){
  P = t(t(X) + fit_dm(X)[,1])
  P = P / rowSums(P)

  mu = colMeans(coda.base::coordinates(P, 'alr'))
  sigma = stats::cov(coda.base::coordinates(P, 'alr'))
  d = length(mu)

  EM_ITER = TRUE
  N_ITER = 0
  while(EM_ITER & N_ITER < 100){
    N_ITER = N_ITER + 1
    res = list()
    for (i in 1:nrow(X)) {
      x = X[i, ]
      m = coda.base::coordinates(P[i,], 'alr')
      V = rep(1, length(m))
      suppressWarnings(
        opt_pars <- optim(c(m,V), F_optim, gr = gr_optim, x, mu, sigma, control = list(maxit = 10000))
      )
      if(opt_pars$convergence != 0){
        warning(sprintf("Not convergence. Error %d", opt_pars$convergence))
      }

      m = opt_pars$par[+(1:d)]
      V = opt_pars$par[-(1:d)]
      res[[i]] = list(mu = m, sigma = V)
    }

    mu_new = rowMeans(sapply(res, function(x) x$mu))
    sigma_new = matrix(rowMeans(sapply(res, function(x) diag(x$sigma) + (x$mu - mu_new) %*% t(x$mu - mu_new))), ncol = d)

    tol = max(abs(sigma - sigma_new))
    EM_ITER = tol > 1e-03
    # print(tol)
    # print(mu)
    mu = mu_new
    sigma = sigma_new
  }

  P.rpl = coda.base::composition(t(sapply(res, function(x)x$mu)), 'alr')
  P.rpl
}
