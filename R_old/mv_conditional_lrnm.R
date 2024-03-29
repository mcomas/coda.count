sinv = function(S){
  chol2inv(chol(S))
}

#' @export
fit_mv_conditional_lrnm = function(X, B1 = coda.base::ilr_basis(ncol(X)), probs = FALSE, sobol.d_length = 100,
                                   eps = 0.001, max_iter = 500){
  Y = X
  isZero = X==0
  Y[isZero] = 0.66

  H1 = coda.base::coordinates(Y, B1)
  pca = stats::prcomp(H1, center = FALSE)

  B1sub = pca$rotation

  H1sub = coda.base::coordinates(Y, as.matrix(B1 %*% B1sub))

  Moments = stats::cov.wt(H1sub)

  B = B1 %*% B1sub
  mu = Moments$center
  sigma = Moments$cov
  order = 10
  pars = list(sigma = 0)
  ITER = 1
  mu_prev = Inf
  sigma_prev = Inf
  while(max(abs(mu_prev - mu)) > eps & ITER < max_iter){
    # print(ITER)
    mu_prev = mu
    sigma_prev = sigma
    pars = mv_conditional_lrnm_EMstep(X, mu, sigma, B, sobol.d_length)
    mu = pars$mu
    sigma = pars$sigma
    inv_sigma = sinv(sigma)
    ITER = ITER + 1
  }

  Moments = mv_conditional_lrnm_Estep(X, mu, sigma, B, sobol.d_length)

  H = sapply(Moments, function(m) m[,1+length(mu)]) |> t()

  Ximp = coda.base::composition(H, B)
  return(list(mu = mu, sigma = sigma, P = Ximp, B1sub = B1sub))
}
#' @export
mv_conditional_lrnm_EMstep = function(X, mu, sigma, B, sobol.d_length = 100){
  d = length(mu)
  moments = mv_conditional_lrnm_Estep(X, mu, sigma, B, sobol.d_length)
  moments_sum = Reduce(`+`, moments)
  mu = moments_sum[,d+1] / nrow(X)
  sigma = moments_sum[,1:d] / (nrow(X)) - mu %*% t(mu)
  list(mu = mu, sigma = sigma)
}

#' @export
mv_conditional_lrnm_Estep = function(X, mu, sigma, B, sobol.d_length = 100){
  zPatterns = apply(X==0, 1, function(i) paste(as.integer(i), collapse=''))
  lPatterns = split(1:nrow(X), zPatterns)
  lMoments = lapply(lPatterns, function(I){
    # print(I)
    # degenerate_conditional_trunc_posterior_moments_hermite_one_pattern
    mv_conditional_lrnm_Estep_OnePattern(X[I,,drop=FALSE], mu, sigma, B, X[I[1],] == 0, sobol.d_length)
    # Xs = X[I,,drop=FALSE]; mu1 = mu; sigma1 = sigma; B1 = B; iZ = X[I[1],] == 0
  })
  lMomentsOrd = list()
  for(i in 1:length(lPatterns)){
    lMomentsOrd[lPatterns[[i]]] = lMoments[[i]]
  }
  lMomentsOrd
}

# posterior_moments_hermite = function(X, mu, sigma, B){
#   d = length(mu)
#   inv_sigma = chol2inv(chol(sigma))
#   Binv = t(MASS::ginv(B))
#   lapply(1:nrow(X), function(i){
#     x = X[i,]
#     N_approx_si = c_posterior_approximation_vec_sigma_inverse(x, mu, inv_sigma, Binv)
#     c_moments_lrnm_hermite_sigma_inverse(x, N_approx_si[,d+1], N_approx_si[,-d-1], mu, inv_sigma, Binv, order = order, mu_centering = rep(0, d))
#   })
# }

mv_conditional_lrnm_Estep_OnePattern = function(Xs, mu1, sigma1, B1, iZ, sobol.d_length){
  sZ = sum(iZ)
  sNZ = ncol(Xs) - sZ
  if(sZ == 0){
    lMom = apply(coda.base::coordinates(Xs, B1), 1, function(h){
      unname(cbind(h %*% t(h), h))
      # M2 = h %*% t(h)
      # eig = eigen(M2)
      # unname(cbind(eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors), h))
    }, simplify = FALSE)
    return(lMom)
  }

  B2 = matrix(0, nrow = ncol(Xs), ncol = ncol(Xs)-1)
  B2[,sZ] = as.matrix(coda.base::sbp_basis(matrix(2*iZ - 1, ncol = 1), silent = TRUE))
  if(sZ>1){
    BZ = as.matrix(coda.base::ilr_basis(sZ))
    B2[iZ,2:sZ-1] = BZ
  }
  if(sNZ>1){
    BnZ = as.matrix(coda.base::ilr_basis(sNZ))
    B2[!iZ,sZ+2:sNZ-1] = BnZ
    H2 = log(Xs[,!iZ]) %*% BnZ
  }

  # Change of basis: B1 -> B2
  Bt = MASS::ginv(B1) %*% B2
  inv_Bt = MASS::ginv(Bt)

  mu2 = mu1 %*% Bt
  sigma2 = t(Bt) %*% sigma1 %*% Bt

  inv_sigma2 = chol2inv(chol(sigma2))
  # inv_sigma2 = MASS::ginv(sigma2)
  inv_B2 = t(MASS::ginv(B2))
  d = ncol(B2)
  if(sNZ == 1){
    Z = randtoolbox::sobol((d+1) * sobol.d_length, dim = d, normal = TRUE) |> t()
    return(lapply(1:nrow(Xs), function(i){
      x = Xs[i,]

      N_approx_si = c_lrnm_posterior_approximation_vec_sigma_inverse(x, mu2, inv_sigma2, inv_B2)
      MomentsH1 = c_moments_lrnm_montecarlo_sigma_inverse(x, N_approx_si[,d+1], N_approx_si[,1:d,drop=FALSE],
                                                          mu2, inv_sigma2,
                                                          inv_B2, Z = Z, mu_centering = rep(0,d))


      M1 = as.vector(MomentsH1[,d+1] %*% inv_Bt)
      M2 = t(inv_Bt) %*% MomentsH1[,1:d,drop=FALSE] %*% inv_Bt
      unname(cbind(M2, M1))
    }))
  }

  i1 = 1:sZ
  i2 = -i1

  inv_sigma2_i2 = sinv(sigma2[i2,i2])
  # inv_sigma2_i2 = MASS::ginv(sigma2[i2,i2])
  # inv_sigma2_i2b = pinv_sympd(sigma2[i2,i2])
  Z = randtoolbox::sobol((sZ+1) * sobol.d_length, dim = sZ, normal = TRUE) |> t()
  lapply(1:nrow(Xs), function(i){
    h2 = H2[i,]
    x = Xs[i,]

    mu_c = as.vector(mu2[i1] + sigma2[i1,i2] %*% inv_sigma2_i2 %*% (h2-mu2[i2]))
    sigma_c = sigma2[i1,i1] - sigma2[i1,i2] %*% inv_sigma2_i2 %*% sigma2[i2,i1]
    inv_sigma_c = sinv(sigma_c)

    N_approx = c_lrnm_cond_posterior_approximation_vec(x, mu_c, inv_sigma_c, h2, inv_B2)
    MomentsH1 = c_moments_lrnm_cond_montecarlo_sigma_inverse(x, N_approx[,i2], sinv(N_approx[,i1,drop=FALSE]),
                                                             mu_c, inv_sigma_c, h2,
                                                             inv_B2, Z, mu_centering = rep(0,sZ))

    # print(as.vector(eigen(rbind(
    #   cbind(MomentsH1[,i1,drop=FALSE], MomentsH1[,i2] %*% t(h2)),
    #   cbind(h2 %*% t(MomentsH1[,i2]), h2 %*% t(h2))))$values))
    M1 = as.vector(c(MomentsH1[,i2], h2) %*% inv_Bt)
    M2 = t(inv_Bt) %*% rbind(
      cbind(MomentsH1[,i1,drop=FALSE], MomentsH1[,i2] %*% t(h2)),
      cbind(h2 %*% t(MomentsH1[,i2]), h2 %*% t(h2))) %*% inv_Bt
    # print(eigen(M2)$values)
    # M2 = Matrix::as.matrix(Matrix::nearPD(M2)$mat)
    unname(cbind(M2, M1))
  })

}

