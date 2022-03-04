library(coda.base)
library(coda.count)
library(zCompositions)
library(mvtnorm)
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_parliament-seed_00002"

###############
load(sprintf("/home/marc/Research/LRNM-REPLACEMENT/sim-01/data/%s.RData", GEN))

d = ncol(X) - 1
M = coordinates(colSums(X))
S = diag(d)

B0 = matrix(0, nrow = 1+length(M), ncol = length(M))
Bd = lapply(1:ncol(X), ilr_basis)


iZ = X == 0
sZ = rowSums(iZ)
sNZ = ncol(X) - sZ

wZ = which(sZ > 0 & sNZ > 1)
wNZ1 = which(sNZ == 1)

lB = lapply(wZ, function(i){
  B = B0
  B[,1] = sbp_basis(matrix(2*iZ[i,] - 1, ncol = 1), silent = TRUE)
  if(sZ[i] > 1){
    B[iZ[i,],2:sZ[i]] = Bd[[sZ[i]]]
  }
  B[!iZ[i,],-(1:sZ[i])] = Bd[[sNZ[i]]]
  return(B)
})

lBt = lapply(lB, function(B) t(B) %*% Bd[[1+d]])

########## Start iteration here
IT = 0
CONT = TRUE
loglikN_prev = NA

IT = IT + 1

lMt = lapply(lBt, function(Bt) Bt %*% M)
lSt = lapply(lBt, function(Bt) Bt %*% S %*% t(Bt))

lh1 = lapply(wZ, function(i) rep(0, sZ[i]))
lh2 = lapply(wZ, function(i) as.vector(t(Bd[[sNZ[i]]]) %*% log(X[i,][!iZ[i,]])))

i = 1
x = X[wZ[i],]
B = lB[[i]]
Mt = lMt[[i]]
St = lSt[[i]]

h1 = lh1[[i]]
h2 = lh2[[i]]

eps = rep(1, sZ[wZ[i]])

i1 = 1:sZ[wZ[[i]]]
i2 = -i1
j1 = 1
j2 = -j1

St2.inv = MASS::ginv(St[j2,j2])
Mc = as.vector(Mt[j1] + St[j1,j2,drop=F] %*% St2.inv %*% (c(h1[j2], h2)-Mt[j2]))
Sc = St[j1,j1] - St[j1,j2] %*% St2.inv %*% St[j2,j1]
invSc = MASS::ginv(Sc)

l_indep_lrnm_cond_join_d1(h1, j1-1, x, Mc, invSc, h2, B)
l_lrnm_cond1_join_d1(c(h1,h2), j1-1, x, Mc, invSc, B)

l_indep_lrnm_cond_join_d2(h1, j1-1, x, Mc, invSc, h2, B)
l_lrnm_cond1_join_d2(c(h1,h2), j1-1, x, Mc, invSc, B)

l_indep_lrnm_cond_join_maximum(h1, j1-1, x, Mc, invSc, h2, B)
h = c(h1,h2)
h[j1] = l_lrnm_cond1_join_maximum(h, j1-1, x, Mc, invSc, B)
mu = h[j1]
sigma = -1/l_lrnm_cond1_join_d2(h, j1-1, x, Mc, invSc, B)

MOMENTS = c_moments_lrnm_cond1_hermite(h, 1, x, mu, sigma, Mt, St, B, 10, 0)

dcond = function(h_){
  h_f = h
  sapply(h_, function(h__){
    h_f[j1] = h__
    dmvnorm(h_f[j2], Mt[j2], St[j2,j2]) * dnorm(h_f[j1], Mc, sqrt(Sc[1,1])) * dmultinom(x, prob = composition(h_f, B))
  })
}
eps = 0.1
hs = seq(-10,1,eps)
dhs = dcond(hs)
dhs = dhs / sum(dhs * eps)
plot(hs, dhs, type = 'l', lwd = 3)
abline(v = mu, col = 'blue')
points(hs, dnorm(hs, mu, sqrt(sigma)), type = 'l', col = 'blue')
abline(v = MOMENTS[2], col = 'red')
points(hs, dnorm(hs, MOMENTS[2], sqrt(MOMENTS[1] - MOMENTS[2]^2)), type='l', col='red')


moments = matrix(0, ncol = 2, nrow = sZ[wZ[i]])
while(max(eps) > 1e-8){

  for(j1 in 1:sZ[wZ[i]]){
    j2 = which(i1 != j1)
    mom = c_moments_indep_lrnm_cond_hermite(X[wZ[i],], j1, napprox[j1,sZ[wZ[i]]+1], napprox[j1,j1], Mt, St, c(h1[j2],h2), B, mu_centering = 0, order = 1000)

    eps[j1] = abs(mom[,2] - h1[j1])
    h1[j1] = mom[,2]
    moments[j1, ] = mom
  }
  # print(moments)
}
