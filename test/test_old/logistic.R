library(coda.base)
logistic_alr = function(h) factorial(length(h)) * prod(composition(h,'alr'))
logistic_ilr = function(h) factorial(length(h)) * sqrt(length(h)+1) * prod(composition(h))

logistic_alr(0)
logistic_ilr(0)


h = seq(-5,5,0.1)
plot(h, sapply(h, logistic_alr), type = 'l')
NLALR = c_logistic_alr_approximation(1)
points(h, dnorm(h, NLALR[2], sqrt(NLALR[1])), type = 'l', col = 'red')


plot(h, sapply(h, logistic_ilr), type = 'l')
NLILR = c_logistic_ilr_approximation(1)
points(h, dnorm(h, NLILR[2], sqrt(NLILR[1])), type = 'l', col = 'red')


library(cubature)
pcubature(logistic_alr, lowerLimit = c(-15), upperLimit = c(15))
pcubature(logistic_ilr, lowerLimit = c(-15), upperLimit = c(15))
# 1
pcubature(logistic, lowerLimit = c(-15,-15), upperLimit = c(15,15))
# 1/2
vegas(logistic, lowerLimit = c(-5,-5,-5), upperLimit = c(5,5,5))
# 1/(2*3)
vegas(logistic, lowerLimit = c(-15,-15,-15,-15), upperLimit = c(15,15,15,15))
# 1/(2*3*4)

logistic_max = function(d){
  ln_Mx = sum(log(2:d)) - (d+0.5) * log(d+1)
  exp(ln_Mx)
}
logistic_max(10)
logistic(rep(0,10))

logistic_sigma = function(d){
  ln_Mx = sum(log(2:d)) - (d+0.5) * log(d+1)
  ln_sigma = 2/d * ln_Mx - log(2*pi)
  exp(ln_sigma)
}
logistic_sigma(1)
