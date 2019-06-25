library(coda.base)
library(coda.count)

library(gtools)

ALPHA = c(10,100)
dirichlet_lebesgue = function(p) ddirichlet(cbind(p,1-p), ALPHA)
integrate(dirichlet_lebesgue, 0, 1)

h = seq(-10, 10, 0.01)
p = composition(cbind(h))[,1]

dirichlet_aitchison = function(h){
  P = composition(cbind(h))
  ddirichlet(P, ALPHA) * sqrt(ncol(P)) * apply(P, 1, prod)
}
integrate(dirichlet_aitchison, -100, 100)

plot(h, dirichlet_lebesgue(p), type = 'l')
points(h, dirichlet_aitchison(h), type = 'l', col='blue')


plot(p, dirichlet_lebesgue(p), type = 'l')
points(p, dirichlet_aitchison(h), type = 'l', col='blue')
