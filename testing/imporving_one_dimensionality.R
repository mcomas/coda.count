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
SIGMA = M %*% diag(rlnorm(d)) %*% t(M)

mu0 = sample(c(-0.5,0,0.5), d, replace = TRUE)
sigma0 = SIGMA
B0 = ilr_basis(d+1)


X = rlrnm(100, 50, mu0, sigma0, B0)

i = which.max(rowSums(X == 0))
x = X[i,]
n0 = sum(x == 0)

B = sbp_basis(2*cbind(x == 0)-1, fill = TRUE)
colnames(B) = sprintf("b%d", 1:d)
rownames(B) = sprintf("x%d", 1:(d+1))
plot_balance(B, main = sprintf("B(x=[%s])", paste(x, collapse=',')))

Bt = MASS::ginv(B0) %*% B
mu_B = mu0 %*% Bt
sigma_B = t(Bt) %*% sigma0 %*% Bt


I1 = 1
I2 = -(1:n0)

h2 = coordinates(x, B[,I2])

inv_sigma_B_I2 = pinv_sympd(sigma_B[I2,I2])
mu_c = as.vector(mu_B[I1] + sigma_B[I1,I2] %*% inv_sigma_B_I2 %*% (h2-mu_B[I2]))
sigma_c = sigma_B[I1,I1] - sigma_B[I1,I2] %*% inv_sigma_B_I2 %*% sigma_B[I2,I1]
inv_sigma_c = pinv_sympd(sigma_c)

d1 = length(I1)
N_h1.x = c_lrnm_cond_posterior_approximation_vec(x, mu_c, inv_sigma_c, h2, cbind(B[,I1], B[,I2]))
mu_h1.x = N_h1.x[,d1+1]
sigma_h1.x = N_h1.x[1:d1,1:d1,drop=FALSE]


Z = matrix(rnorm(1000*d1), nrow = d1, ncol = 1000)
moments_h1.x = c_moments_lrnm_cond_montecarlo_sigma_inverse(x, mu_h1.x, pinv_sympd(sigma_h1.x),
                                                            mu_c, inv_sigma_c, h2,
                                                            cbind(B[,I1], B[,I2]), Z, mu_centering = rep(0,d1))
moments_h1.x[,d1+1]
moments_h1.x[1:d1,1:d1,drop=FALSE]

x/composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))
composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))[x==0]

# When I1 = 1:n0
#
# > x/composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))
#       x1       x2       x3       x4       x5       x6       x7       x8       x9      x10      x11
# 51.44938  0.00000  0.00000  0.00000  0.00000 51.44938  0.00000  0.00000  0.00000  0.00000 51.44938
# > composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))[x==0]
#           x2           x3           x4           x5           x7           x8           x9          x10
# 0.0041431095 0.0034094285 0.0055928956 0.0015893392 0.0036707773 0.0008601778 0.0027059679 0.0061992328
#
#
# When I1 = 1
#
# > x/composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))
#       x1       x2       x3       x4       x5       x6       x7       x8       x9      x10      x11
# 51.57835  0.00000  0.00000  0.00000  0.00000 51.57835  0.00000  0.00000  0.00000  0.00000 51.57835
# > composition(c(moments_h1.x[,d1+1], h2), cbind(B[,I1], B[,I2]))[x==0]
#          x2          x3          x4          x5          x7          x8          x9         x10
# 0.003825121 0.003825121 0.003825121 0.003825121 0.003825121 0.003825121 0.003825121 0.003825121
