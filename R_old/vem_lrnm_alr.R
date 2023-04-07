optimise_xi_m_V_alr = function(m, V, xi, x, mu, sigma){
  xi = xi_optimise_alr(m, V, xi, x, mu, sigma)

  F_alr_m = function(m){
    -F_alr(m, V, xi, x, mu, sigma)
  }
  gr_m = function(m){
    DF_m_alr(m, V, xi, x, mu, sigma)
  }
  m = optim(m, F_alr_m, gr_m)$par

  F_alr_V = function(V){
    -F_alr(m, V^2, xi, x, mu, sigma)
  }
  gr_V = function(V){
    d1v_alr(m, V, xi, x, mu, sigma)
  }
  V = (optim(sqrt(V), F_alr_V, gr_V)$par)^2

  list(xi = xi, m = m, V = V)
}
m_Newton_Step_alr = function(m, V, xi, x, mu, sigma){
  Der = DF_m_alr(m, V, xi, x, mu, sigma)
  Hes = HF_m_alr(m, V, xi, x, mu, sigma)
  m - 0.1 * solve(Hes) %*% Der
}
V_Newton_Step_alr = function(m, V, xi, x, mu, sigma){
  d1 = d1v_alr(m, V, xi, x, mu, sigma)
  d2 = d2v_alr(m, V, xi, x, mu, sigma)
  V - 0.1 * d1/d2
}
#' @export
fit_vem_lrnm_alr = function(X){
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
      xi = 1
      EPS = 1
      ITER = 1
      while (ITER < 100 & EPS > 1e-05) {
        opt_pars = optimise_xi_m_V_alr(m, V, xi, x, mu, sigma)
        # print(opt_pars$m)
        ITER = ITER + 1
        EPS = max(abs(m - opt_pars$m))
        xi = opt_pars$xi
        m = opt_pars$m
        V = opt_pars$V
      }
      res[[i]] = list(mu = m, sigma = V, iter = ITER)
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
