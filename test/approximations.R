library(coda.base)
library(coda.count)
dirichlet_gaussian_approx = function(ALPHA){
  D = length(ALPHA)
  ILR2ALR = t(ilr_basis(D)) %*% alr_basis(D)
  ALR2ILR = solve(ILR2ALR)

  MU_ALR = digamma(ALPHA[-D]) - digamma(ALPHA[D])
  MU = coordinates(composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR

  SIGMA_ALR = diag(trigamma(ALPHA[-D]), nrow = D-1) + trigamma(ALPHA[D])
  SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR
  return(list(MU = MU, SIGMA = SIGMA))
}
gaussian_product = function(MU1, SIGMA1, MU2, SIGMA2){
  invSIGMA12 = solve(SIGMA1+SIGMA2)
  SIGMA3 = SIGMA1 %*% invSIGMA12 %*% SIGMA2
  MU3 = SIGMA2 %*% invSIGMA12 %*% t(t(MU1)) + SIGMA1 %*% invSIGMA12 %*% t(t(MU2))
  list(MU = MU3, SIGMA = SIGMA3)
}
gaussian_division = function(MU2, SIGMA2, MU3, SIGMA3){
  invSIGMA3 = solve(SIGMA3)
  invSIGMA2 = solve(SIGMA2)
  SIGMA1 = solve(invSIGMA3 - invSIGMA2)
  MU1 = SIGMA1 %*% invSIGMA3 %*% MU3 - SIGMA1 %*% invSIGMA2 %*% MU2
  list(MU = MU1, SIGMA = SIGMA1)
}
