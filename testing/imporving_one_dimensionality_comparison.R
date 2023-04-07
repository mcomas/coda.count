library(coda.count)
library(coda.base)

set.seed(1)
d = 10

rotation = random_rotation_matrix_incl_flip <- function(n){
  QR <- qr(matrix(rnorm(n^2), ncol=n))          # A = QR
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))  # ensure diag(R) > 0
  return(M)
}
M = rotation(d)
SIGMA = M %*% diag(0.1*rlnorm(d)) %*% t(M)

mu0 = sample(c(-0.5,0,0.5), d, replace = TRUE)
sigma0 = SIGMA
B0 = ilr_basis(d+1)

X = rlrnm(1000, 50, mu0, sigma0, B0)

mu_ = mu0 + 1
mu_next = mu0
sigma_next = sigma0

i0 = rowSums(X == 0) > 0
zPatterns = apply(1L*(X == 0), 1, paste, collapse = '')
zPatterns_unique = unique(zPatterns[i0])

H = coordinates(X[!i0,,drop=FALSE], B0)
sumMoments_nonzero = cbind(t(H) %*% H, colSums(H))

if(sum(i0) > 0){
  while( max(abs(mu_next-mu_)) > 0.001){

    mu_ = mu_next
    sigma_ = sigma_next

    lMoments = lapply(zPatterns_unique, function(zPattern){
      X_ = X[zPatterns == zPattern,,drop=FALSE]

      x = X_[1,]
      n0 = sum(x == 0)

      B = sbp_basis(2*cbind(x == 0)-1, fill = TRUE)
      if(FALSE){
        colnames(B) = sprintf("b%d", 1:d)
        rownames(B) = sprintf("x%d", 1:(d+1))
        plot_balance(B, main = sprintf("B(x=[%s])", paste(x, collapse=',')))
      }
      Bt = MASS::ginv(B0) %*% B


      mu_B = mu_ %*% Bt
      sigma_B = t(Bt) %*% sigma_ %*% Bt

      I1 = 1:n0
      I2 = -(1:n0)

      H2 = coordinates(X_, B[,I2])

      inv_sigma_B_I2 = pinv_sympd(sigma_B[I2,I2])
      sigma_c = sigma_B[I1,I1] - sigma_B[I1,I2] %*% inv_sigma_B_I2 %*% sigma_B[I2,I1]

      eig = eigen(sigma_c)
      h_pc1 = B[,I1,drop=FALSE] %*% eig$vectors[,1]
      # h_pc1
      B_ = cbind(h_pc1, B[,I2])
      inv_B0_ = MASS::ginv(MASS::ginv(B0) %*% B_)
      # t(B_) %*% B_


      MU_c = apply(H2, 1, function(h2) as.vector(mu_B[I1] + sigma_B[I1,I2] %*% inv_sigma_B_I2 %*% (h2-mu_B[I2])))

      MU_pc1 = as.vector(eig$vectors[,1] * MU_c)
      sigma_pc1 = t(eig$vectors[,1,drop=FALSE]) %*% sigma_c %*% eig$vectors[,1,drop=FALSE]
      inv_sigma_pc1 = pinv_sympd(sigma_pc1)

      d1 = 1
      # Z = matrix(rnorm(1000*d1), nrow = d1, ncol = 1000)
      Moments = lapply(1:nrow(X_), function(i){
        N_h1.x = c_lrnm_cond_posterior_approximation_vec(X_[i,], MU_pc1[i], inv_sigma_pc1, H2[i,], cbind(h_pc1, B[,I2]))
        # MomentsH1 = c_moments_lrnm_cond_montecarlo_sigma_inverse(X_[i,], N_h1.x[,d1+1], pinv_sympd(N_h1.x[1:d1,1:d1,drop=FALSE]),
        #                                                          MU_pc1[i], inv_sigma_pc1, H2[i,],
        #                                                          cbind(h_pc1, B[,I2]), Z, mu_centering = rep(0,d1))
        MomentsH1 = c_moments_lrnm_cond_hermite_1d(X_[i,], N_h1.x[,d1+1], N_h1.x[1:d1,1:d1,drop=FALSE],
                                                   MU_pc1[i], inv_sigma_pc1, H2[i,],
                                                   cbind(h_pc1, B[,I2]), 100, mu_centering = rep(0,d1))
        list(
          M1 = as.vector(c(MomentsH1[,2], H2[i,]) %*% inv_B0_),
          M2 = t(inv_B0_) %*% rbind(
            cbind(MomentsH1[,1], MomentsH1[,2] %*% t(H2[i,])),
            cbind(H2[i,] %*% t(MomentsH1[,2]), H2[i,] %*% t(H2[i,]))) %*% inv_B0_)

      })


      cbind(Reduce(`+`, lapply(Moments, function(m) m$M2)),
            Reduce(`+`, lapply(Moments, function(m) m$M1)))

    })
    sumMoments = sumMoments_nonzero + Reduce(`+`, lMoments)
    mu_next = sumMoments[,d+1] / nrow(X)
    sigma_next = sumMoments[,1:d, drop=FALSE] / nrow(X) - mu_next %*% t(mu_next)
  }
}else{
  sumMoments = sumMoments_nonzero
  mu_next = sumMoments[,d+1] / nrow(X)
  sigma_next = sumMoments[,1:d, drop=FALSE] / nrow(X) - mu_next %*% t(mu_next)
}
plot(mu0, ylim = c(-1,1))
points(mu_next, col = 'blue')

