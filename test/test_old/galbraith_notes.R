library(coda.count)

set.seed(1)

n = 100
MU = -0.5
SIGMA = 1.5

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

log(nu/(1-nu))
