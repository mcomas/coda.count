optimise_xi_m_V = function(m, V, xi, x, mu, sigma, B){
  xi = xi_optimise(m, V, xi, x, mu, sigma, B)


  m = m_Newton_Step(m, V, xi, x, mu, sigma, B)



  V = V_Newton_Step(m, V, xi, x, mu, sigma, B)
  list(xi = xi, m = m, V = V)
}
m_Newton_Step = function(m, V, xi, x, mu, sigma, B){
  Der = DF_m(m, V, xi, x, mu, sigma, B)
  Hes = HF_m(m, V, xi, x, mu, sigma, B)
  m - solve(Hes) %*% Der
}
V_Newton_Step = function(m, V, xi, x, mu, sigma, B){
  d1 = d1v(m, V, xi, x, mu, sigma, B)
  d2 = d2v(m, V, xi, x, mu, sigma, B)
  V - 0.1 * d1/d2
}
#' @export
fit_vem_lrnm = function(X, B1 = coda.base::ilr_basis(ncol(X)), probs = FALSE, hermite.order = 5,
                        eps = NULL, max_iter = 500){
  P = t(t(X) + fit_dm(X)[,1])
  P = P / rowSums(P)
  B = coda.base::alr_basis(ncol(X))

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
        opt_pars = optimise_xi_m_V(m, V, xi, x, mu, sigma, B)
        # print(opt_pars$V)
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
}
