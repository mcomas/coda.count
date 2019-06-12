ALPHA = c(1,3,2)

MU_ALR = digamma(ALPHA[1:2])-digamma(ALPHA[3])
SIGMA_ALR = diag(trigamma(ALPHA[1:2])) + trigamma(ALPHA[3])

inv.trigamma = function(y){
  x <- 0.5+1/y
  tri <- trigamma(x)
  x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
  tri <- trigamma(x)
  x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
  tri <- trigamma(x)
  x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
  tri <- trigamma(x)
  x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
  x
}
M = c(1,2)
M_COMP = composition(M)



S = diag(2) * 3
S_ALR = (t(solve(ALR2ILR)) %*% S %*% solve(ALR2ILR))
inv.trigamma(S_ALR[1,1] / 2 )


ALR2ILR = MASS::ginv(alr_basis(3)) %*% ilr_basis(3)
(MU = coordinates(composition(MU_ALR, 'alr')))
(SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR)

ID12 = t(solve(ALR2ILR)) %*% (diag(2)*3) %*% solve(ALR2ILR)
ID12
t(ALR2ILR) %*% ID12 %*% ALR2ILR

t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR
