
#####
#' @export
fit_lrnm_vem_original = function(X){
  F_original = function(xi, m, V, x, mu_, sigma_, B){
    sigma_inv = pinv_sympd(sigma_)
    S1 = + sum(x * m)
    S2 = - sum(x) * ((1/xi) * sum( exp(m+V/2) ) - 1 + log(xi))
    S3 = - 0.5 * log( det(t(B) %*% sigma_ %*% B) )
    S4 = - 0.5 * sum(((m - mu_) %*% sigma_inv) * (m - mu_) )
    S5 = - 0.5 * sum(diag(sigma_inv) * V)
    S6 = + 0.5 * sum(V) + 0.5 * d
    S1 + S2 + S3 + S4 + S5 + S6
  }
  gr_m_F_original = function(xi, m, V, x, mu_, sigma_, B){
    sigma_inv = pinv_sympd(sigma_)
    as.vector(x - sigma_inv %*% (m-mu_) - sum(x) / xi * exp(m+0.5*V))[1:d]
  }
  gr_v_F_original = function(xi, m, V1, x, mu_, sigma_, B){
    sigma_inv = pinv_sympd(sigma_)
    (1/V1 - V1 * diag(sigma_inv) - sum(x) * (1/xi) * exp(m + 0.5 * V1^2) * V1)[1:d]
  }

  D = ncol(X)
  d = D - 1
  B = rbind(diag(d), 0)

  mu_new = rep(0, D)
  sigma_new = cbind(rbind(diag(d), 0), 0)

  mu_ = c(rep(1,d), 0)

  while(max(abs(mu_ - mu_new)) > 0.0001){
    mu_ = mu_new
    sigma_ = sigma_new

    pars = lapply(1:nrow(X), function(x) list(m = c(rep(1,d), 0),
                                              V = c(rep(1,d), 0)))


    for(I in 1:nrow(X)){

      xi = sum(exp(pars[[I]]$m+0.5*pars[[I]]$V))
      # F_original(xi, m, V, x, mu_, sigma_, B)

      m_opt = optim(pars[[I]]$m[1:d],
                    fn = function(m_) -F_original(xi, c(m_, 0), pars[[I]]$V, X[I,], mu_, sigma_, B),
                    gr = function(m_) -gr_m_F_original(xi, c(m_, 0), pars[[I]]$V, X[I,], mu_, sigma_, B))

      pars[[I]]$m[1:d] = m_opt$par
      # F_original(xi, m, V, x, mu_, sigma_, B)

      V_opt = optim(sqrt(pars[[I]]$V[1:d]),
                    fn = function(V_) -F_original(xi, pars[[I]]$m, c(V_,0)^2, X[I,], mu_, sigma_, B),
                    gr = function(V_) -gr_v_F_original(xi, pars[[I]]$m, c(V_, 0), X[I,], mu_, sigma_, B), method = 'L-BFGS-B')
      pars[[I]]$V[1:d] = V_opt$par^2
    }
    mu_new = rowMeans(sapply(pars, function(p) p$m))
    sigma_new = matrix(rowMeans(sapply(pars, function(p) diag(p$V) + (p$m-mu_) %*% t(p$m-mu_))),
                       nrow = D, ncol = D)
  }
  list(m=mu_new[1:d],
       sigma=sigma_new[1:d,1:d],
       P=composition(t(sapply(pars, function(p) p$m))[,1:d], 'alr'))
}
# F_original(xi, m, V, x, mu_, sigma_, B)
