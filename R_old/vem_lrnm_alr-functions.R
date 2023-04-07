F_alr = function(m, V, xi, x, mu, sigma){
  invSigma = MASS::ginv(sigma)
  as.vector(sum(x[-length(x)] * m) - sum(x) * (xi^-1 * sum(exp(m + 0.5 * V))- 1 + log(xi) ) +
      - 0.5 * (log(det(sigma)) + t(m-mu) %*% invSigma %*% (m-mu) + sum(diag(invSigma) * V)) +
      + 0.5 * (sum(log(V)) + length(m)))
}

DF_m_alr = function(m, V, xi, x, mu, sigma){
  invSigma = MASS::ginv(sigma)
  x[-length(x)] - as.vector(invSigma %*% (m-mu)) -
    as.vector(sum(x) * xi^-1 * exp(m + 0.5 * V))
}

HF_m_alr = function(m, V, xi, x, mu, sigma){
  invSigma = MASS::ginv(sigma)
  E = diag(exp(m + 0.5 * V))
  -invSigma - sum(x) * xi^-1 * E
}

d1v_alr = function(m, V, xi, x, mu, sigma){
  invSigma = MASS::ginv(sigma)
  # colSums(Binv * Binv * as.vector(E)) * V == (t(Binv * Binv) * V) %*% E
  as.vector(1/V - V * diag(invSigma) - sum(x) * xi^-1 * exp(m + 0.5 * V^2) * V)
}

d2v_alr = function(m, V, xi, x, mu, sigma){
  invSigma = MASS::ginv(sigma)
  E = exp(m + 0.5 * V)
  as.vector(-1/V - diag(invSigma) - sum(x) * xi^-1 * (1 + V) * E )
}

xi_optimise_alr = function(m, V, xi, x, mu, sigma){
  sum(exp(m + 0.5 * V))
}

