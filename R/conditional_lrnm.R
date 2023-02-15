#' @export
fit_conditional_lrnm = function(X, B1 = coda.base::ilr_basis(ncol(X)), probs = FALSE, hermite.order = 5,
                                eps = 1e-4, max_iter = 500){
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
  while(max(abs(mu_prev - mu)) > eps){
    print(ITER)
    mu_prev = mu
    sigma_prev = sigma
    pars = conditional_lrnm_EMstep(X, mu, sigma, B, hermite.order)
    mu = pars$mu
    sigma = pars$sigma
    ITER = ITER + 1
  }

  Moments = conditional_lrnm_Estep(X, mu, sigma, B, hermite.order)

  H = sapply(Moments, function(m) m[,1+length(mu)]) |> t()

  Ximp = coda.base::composition(H, B)
  return(list(mu = mu, sigma = sigma, P = Ximp, B1sub = B1sub))
}
#' @export
conditional_lrnm_EMstep = function(X, mu, sigma, B, hermite.order = 10){
  d = length(mu)
  moments = conditional_lrnm_Estep(X, mu, sigma, B, hermite.order)
  moments_mean = Reduce(`+`, moments) / nrow(X)
  mu = moments_mean[,d+1]
  sigma = moments_mean[,1:d] - mu %*% t(mu)
  list(mu = mu, sigma = sigma)
}

#' @export
conditional_lrnm_Estep = function(X, mu, sigma, B, hermite.order = 10){
  zPatterns = apply(X==0, 1, function(i) paste(as.integer(i), collapse=''))
  lPatterns = split(1:nrow(X), zPatterns)
  lMoments = lapply(lPatterns, function(I){
    # degenerate_conditional_trunc_posterior_moments_hermite_one_pattern
    conditional_lrnm_Estep_OnePattern(X[I,,drop=FALSE], mu, sigma, B, X[I[1],] == 0, hermite.order)
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

conditional_lrnm_Estep_OnePattern = function(Xs, mu1, sigma1, B1, iZ, hermite.order){
  sZ = sum(iZ)
  sNZ = ncol(Xs) - sZ
  if(sZ == 0){
    lMom = apply(coda.base::coordinates(Xs, B1), 1, function(h){
      unname(cbind(h %*% t(h), h))
    }, simplify = FALSE)
    return(lMom)
  }

  B2 = matrix(0, nrow = ncol(Xs), ncol = sNZ)
  B2[,1] = as.matrix(coda.base::sbp_basis(matrix(2*iZ - 1, ncol = 1), silent = TRUE))
  if(sNZ>1){
    BnZ = as.matrix(coda.base::ilr_basis(sNZ))
    B2[!iZ,2:sNZ] = BnZ
    H2 = log(Xs[,!iZ]) %*% BnZ
  }

  # Change of basis: B1 -> B2
  Bt = MASS::ginv(B1) %*% B2
  inv_Bt = MASS::ginv(Bt)

  mu2 = mu1 %*% Bt
  sigma2 = t(Bt) %*% sigma1 %*% Bt

  # inv_sigma2 = chol2inv(chol(sigma2))
  inv_sigma2 = MASS::ginv(sigma2)
  inv_B2 = t(MASS::ginv(B2))

  if(sNZ == 1){
    return(lapply(1:nrow(Xs), function(i){
      x = Xs[i,]

      N_approx = c_lrnm_posterior_approximation_vec(x, mu2, inv_sigma2, inv_B2)
      MomentsH1 = c_moments_lrnm_hermite(x, N_approx[,2], N_approx[,1,drop=FALSE],
                                         mu2, inv_sigma2,
                                         inv_B2, mu_centering = rep(0,1), order = hermite.order)


      M1 = as.vector(MomentsH1[,2] %*% inv_Bt)
      M2 = t(inv_Bt) %*% MomentsH1[,1] %*% inv_Bt
      unname(cbind(M2, M1))
    }))
  }

  i1 = 1
  i2 = -i1

  # inv_sigma2_i2 = chol2inv(chol(sigma2[i2,i2]))
  inv_sigma2_i2 = MASS::ginv(sigma2[i2,i2])
  # inv_sigma2_i2b = pinv_sympd(sigma2[i2,i2])

  lapply(1:nrow(Xs), function(i){
    h2 = H2[i,]
    x = Xs[i,]

    mu_c = as.vector(mu2[i1] + sigma2[i1,i2] %*% inv_sigma2_i2 %*% (h2-mu2[i2]))
    sigma_c = 1e-10 + abs(sigma2[i1,i1] - sigma2[i1,i2] %*% inv_sigma2_i2 %*% sigma2[i2,i1])
    inv_sigma_c = MASS::ginv(sigma_c)

    N_approx = c_lrnm_cond_posterior_approximation_vec(x, mu_c, inv_sigma_c, h2, inv_B2)
    MomentsH1 = c_moments_lrnm_cond_hermite_1d(x, N_approx[,2], N_approx[,1,drop=FALSE],
                                               mu_c, inv_sigma_c, h2,
                                               inv_B2, mu_centering = rep(0,1), order = hermite.order)


    M1 = as.vector(c(MomentsH1[,2], h2) %*% inv_Bt)
    M2 = t(inv_Bt) %*% rbind(
      cbind(MomentsH1[,1], MomentsH1[,2] %*% t(h2)),
      cbind(h2 %*% t(MomentsH1[,2]), h2 %*% t(h2))) %*% inv_Bt
    unname(cbind(M2, M1))
  })

}
