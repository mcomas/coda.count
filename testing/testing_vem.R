K = 5

x = rmultinom(1, 4, rep(1,K+1))

m = c(rep(1, K), 0)
v2 = rep(0.8,K)
V = diag(c(v2,0))

mu = rep(0.5, K)
mu_ = c(mu, 0)
sigma = toeplitz(K:1/K)
sigma_ = rbind(cbind(sigma,0),0)

xi = 100

B = rbind(diag(K),0)
F_ = function(m, V, xi, x, mu, sigma, B){
  Binv = t(MASS::ginv(B))
  invSigma = MASS::ginv(sigma)
  (x %*% Binv %*% m - sum(x) * (xi^-1 * sum(exp((Binv %*% m) + 0.5 * diag(Binv %*% diag(V) %*% t(Binv)))) - 1 + log(xi) ) +
      - 0.5 * (log(det(sigma)) + t(m-mu) %*% invSigma %*% (m-mu) + sum(diag(invSigma) * V)) +
      + 0.5 * (sum(log(V)) + length(m))) |> as.vector()
}

Rcpp::sourceCpp("src/vem_lrnm.cpp")
c_F(x,m,V,mu_,sigma_,xi,B)
F_(m[1:K],diag(V)[1:K],xi,as.vector(x),mu[1:K],sigma[1:K,1:K],rbind(diag(K),0))

set.seed(6)
X = rlrnm(10, 15, mu, sigma)
str(fit_ <- c_vem_lrnm_fit(t(X), 0.001, 1000))
cbind(t(fit_$m), X)
