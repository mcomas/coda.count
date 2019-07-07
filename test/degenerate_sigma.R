library(coda.count)
library(coda.base)
library(vegan)
data(dune)

Y1 = as.matrix(dune)
alpha1 = fit_dm(Y1)[,1]
H1 = coordinates(t(t(Y1) + alpha1))
X1 = cbind(rep(1, nrow(Y1)))
beta1 = solve(t(X1) %*% X1) %*% t(X1) %*% H1
R1 = H1-X1 %*% beta1
SIGMA1 = (t(R1) %*% R1) / (nrow(X1)-1)
eigen(SIGMA1)$values
barplot(eigen(SIGMA1)$values, col=c(rep(1, nrow(Y1)), rep(2, ncol(Y1)-nrow(Y1))))

Y2 = rmultinomial(100, 10000, c(0.3, 0.2, 0.5))
alpha2 = fit_dm(Y2)[,1]
H2 = coordinates(t(t(Y2) + alpha2), ilr_basis(3, type = 'cdp'))
X2 = cbind(rep(1, nrow(Y2)))
beta2 = solve(t(X2) %*% X2) %*% t(X2) %*% H2
R2 = H2-X2 %*% beta2
SIGMA2 = (t(R2) %*% R2) / (nrow(X2)-1)
eigen(SIGMA2)$value
barplot(eigen(SIGMA2)$values)


SIGMA2
t(eigen(SIGMA2)$vectors) %*% diag(eigen(SIGMA2)$values) %*% eigen(SIGMA2)$vectors

solve(SIGMA2)
eigen(SIGMA2)$vectors %*% diag(1/eigen(SIGMA2)$values) %*% solve(eigen(SIGMA2)$vectors)




Y3 = HardyWeinberg::HWData(100, 100)
alpha3 = fit_dm(Y3)[,1]
H3 = coordinates(t(t(Y3) + alpha3), ilr_basis(3))
X3 = cbind(rep(1, nrow(Y3)))
beta3 = solve(t(X3) %*% X3) %*% t(X3) %*% H
R3 = H3-X3 %*% beta3
SIGMA3 = (t(R3) %*% R3) / (nrow(X3)-1)
eigen(SIGMA3)$value
barplot(eigen(SIGMA3)$values)
