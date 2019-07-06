
library(coda.base)
library(coda.count)
library(ggplot2)

simulation = function(n = 100, MU = -0.5, SIGMA = 1.5){

  beta = rnorm(n, MU, SIGMA)
  theta = exp(beta)/(1+exp(beta))

  m = rep(20, n)
  Ns = rbinom(n, m, theta)
  y = Ns/m

  ######## 1
  z = log(0.5 + Ns) - log(0.5 + m-Ns)

  ######## 2
  sigma = 0.6 * sd(z)
  nu = sum(Ns) / sum(m)

  for(i in 1:20){
    ####### 3
    w = m / (nu * (1-nu) + (m-1) * nu^2 * (1-nu)^2 * sigma^2)

    ####### 4
    sigma = sigma * (sum(w^2*(y-nu)^2) / sum(w) )^0.5
    nu = sum(w * y) / sum(w)
  }

  ## Mu and Sigma estimates according to (Galbraith & Laslett, 1993)
  MU_est1 = log(nu/(1-nu))
  SIGMA_est1 = sigma
  ## Mu estimate according to personal note by Galbraith, 2018
  MU_est2 = log(nu/(1-nu)) + 1/2 * sigma^2 * (nu^2 - (1-nu)^2)


  ########
  X = cbind(Ns, m - Ns)
  mu0 = log(sum(Ns) / (sum(m)-sum(Ns)))
  fit = coda.count::c_lrnm_fit_hermite(X, mu0, matrix(1))
  ILR2ALR = t(ilr_basis(2)) %*% alr_basis(2)
  ## Mu and sigma estimates using hermite integration
  MU_est3 = fit[[1]] %*% ILR2ALR
  SIGMA_est3 = sqrt(t(ILR2ALR) %*% fit[[2]] %*% ILR2ALR)

  # Estimated theta's can be obtained with
  # head(fit[[3]])

  data.frame(
    method = c('Original', 'Galbraith & Laslett, 1993', 'Galbraith, 2018', 'Hermite'),
    mean = c(MU, MU_est1, MU_est2, MU_est3),
    sigma = c(SIGMA, SIGMA_est1, NA_real_, SIGMA_est3)
  )

}

set.seed(1)
REP = replicate(10, simulation(n = 100, MU = -0.5, SIGMA = 1.5), simplify = FALSE)
df = do.call(`rbind`, REP)

boxplot(mean~method, data=df, horizontal=TRUE, las=1, main = 'Mean')
boxplot(sigma~method, data=df, horizontal=TRUE, las=1, main = 'Sigma')


REP = replicate(10, simulation(n = 100, MU = -2, SIGMA = 3), simplify = FALSE)
df = do.call(`rbind`, REP)

boxplot(mean~method, data=df, horizontal=TRUE, las=1, main = 'Mean')
boxplot(sigma~method, data=df, horizontal=TRUE, las=1, main = 'Sigma')
