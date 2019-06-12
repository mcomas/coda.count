library(coda.count)

MU = c(0,0)
SIGMA = diag(2) * log(50)
x = c(1,2,30)

ord = seq(20, 100, 1)
f = sapply(ord, function(order_){
  c_dlrnm_hermite(x = x, mu = MU, sigma = SIGMA, order = order_)
})
plot(ord, f, type = 'l')
