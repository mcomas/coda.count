optimise_xi_m_V = function(m, V, xi, x, mu, sigma, B){
  xi = xi_optimise(m, V, xi, x, mu, sigma, B)

  m_Newton_Step = function(m, V, xi, x, mu, sigma, B){
    Der = DF_m(m, V, xi, x, mu, sigma, B)
    Hes = HF_m(m, V, xi, x, mu, sigma, B)
    m - solve(Hes) %*% Der
  }
  m = m_Newton_Step(m, V, xi, x, mu, sigma, B)

  V_Newton_Step = function(m, V, xi, x, mu, sigma, B){
    d1 = d1v(m, V, xi, x, mu, sigma, B)
    d2 = d2v(m, V, xi, x, mu, sigma, B)
    V - d1/d2
  }

  V = V_Newton_Step(m, V, xi, x, mu, sigma, B)
  list(xi = xi, m = m, V = V)
}

F_ = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  (x %*% Binv %*% m - sum(x) * (xi^-1 * sum(exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv)))) - 1 + log(xi) ) +
      - 0.5 * (log(det(sigma)) + t(m-mu) %*% invSigma %*% (m-mu) + sum(diag(invSigma) * V)) +
      + 0.5 * (sum(log(V)) + length(m))) |> as.vector()
}

DF_m = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  as.vector(x %*% Binv) - as.vector(invSigma %*% (m-mu)) -
    as.vector(sum(x) * xi^-1 * t(Binv) %*% exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv))))
}

HF_m = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  E = exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv)))
  -invSigma - sum(x) * xi^-1 * (t(Binv) %*% (Binv * as.vector(E)))
}

d1v = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  E = exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv)))
  # colSums(Binv * Binv * as.vector(E)) * V == (t(Binv * Binv) * V) %*% E
  1/sqrt(V) - sqrt(V) * diag(invSigma) - sum(x) * xi^-1 * (t(Binv * Binv) * sqrt(V)) %*% E |> as.vector()
}

d2v = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  E = exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv)))
  -1/V - diag(invSigma) - sum(x) * xi^-1 * colSums( t(t(1 + Binv * Binv) * V) * (Binv * Binv * as.vector(E)) )  |> as.vector()
}

xi_optimise = function(m, V, xi_opt, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  sum(exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv))))
}

