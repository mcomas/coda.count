library(coda.count)
library(coda.base)

x = c(0,2,5,0)
M = c(2,1,1)
S = diag(3)

iZ = x==0
sZ = sum(iZ)
sNZ = sum(!iZ)


B = matrix(0, nrow = 1+length(M), ncol = length(M))
B[,1] = sbp_basis(matrix(2*iZ - 1, ncol = 1), silent = TRUE)
if(sZ > 1){
  B[x==0,2:sZ] = ilr_basis(sZ)
}
B[!iZ,-(1:sZ)] = ilr_basis(sNZ)

Bt = t(B) %*% ilr_basis(length(x))

Mt = Bt %*% M
St = Bt %*% S %*% t(Bt)

h2 = t(ilr_basis(sNZ)) %*% log(x[!iZ]) |> as.vector()

i1 = 1:sZ
i2 = -i1

invSt2 = MASS::ginv(St[i2,i2])
Mc = as.vector(Mt[i1] + St[i1,i2] %*% invSt2 %*% (h2-Mt[i2]))
Sc = St[i1,i1] - St[i1,i2] %*% invSt2 %*% St[i2,i1]


napprox = c_lrnm_cond_posterior_approximation_vec(x, Mc, MASS::ginv(Sc), h2, B)

c_moments_lrnm_cond_hermite(x, napprox[,3], napprox[1:2,1:2], Mt, St, h2, B, mu_centering = rep(0,2), order = 5)

j1 = 2
j2 = which(i1 != j1)
c_moments_indep_lrnm_cond_hermite(x, j1, napprox[j1,3], napprox[j1,j1], Mt, St, c(napprox[j2,sZ+1],h2), B, mu_centering = 0, order = 10)


