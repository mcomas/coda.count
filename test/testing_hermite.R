library(coda.count)
library(coda.base)
library(MASS)
D = 4

MU = c(4,3,3)
SIGMA = diag(D-1) * log(50)
SIGMA[1,2] = SIGMA[2,1] = 0.8

B1 = ilr_basis(D)

conv_pars = function(MU, SIGMA, B_to, B_from = ilr_basis(length(MU)+1)){
  A = t(B_from) %*% B_to
  list('MU' = (MU %*% A)[1,],
       'SIGMA' = t(A) %*% SIGMA %*% A)
}

B2 = ilr_basis(D)[D:1,(D-1):1]
p2 = conv_pars(MU, SIGMA, B2, B1)
B3 = alr_basis(D)
p3 = conv_pars(MU, SIGMA, B3, B1)

x = c(1,2,30,2)

ord = seq(2, 100, 1)
f1 = sapply(ord, function(order_){
  c_d_lrnm_hermite(x = x, mu = MU, sigma = SIGMA, t(ginv(B1)), order = order_)
})
f2 = sapply(ord, function(order_){
  c_d_lrnm_hermite(x = x, mu = p2$MU, sigma = p2$SIGMA, t(ginv(B2)), order = order_)
})
f3 = sapply(ord, function(order_){
  c_d_lrnm_hermite(x = x, mu = p3$MU, sigma = p3$SIGMA, t(ginv(B3)), order = order_)
})

ylim_ = range(c(f1,f2,f3))
plot(ord, f1, type = 'l', ylim = ylim_)
points(ord, f2, type='l', col = 'red', lty = 2)
points(ord, f3, type='l', col = 'blue', lty = 3)


ord = seq(2, 100, 2)
e1 = sapply(ord, function(order_){
  c_m1_lrnm_hermite(x = x, mu = MU, sigma = SIGMA, t(ginv(B1)), order = order_)/f1[length(f1)]
})
e2 = sapply(ord, function(order_){
  c_m1_lrnm_hermite(x = x, mu = p2$MU, sigma = p2$SIGMA, t(ginv(B2)), order = order_)/f2[length(f2)]
})
e3 = sapply(ord, function(order_){
  c_m1_lrnm_hermite(x = x, mu = p3$MU, sigma = p3$SIGMA, t(ginv(B3)), order = order_)/f3[length(f3)]
})

tail(t(e1), n = 2)
coordinates(composition(tail(t(e2), n = 2), B2))
coordinates(composition(tail(t(e3), n = 2), 'alr'))

c1 = composition(t(e1), B1)
c2 = composition(t(e2), B2)
c3 = composition(t(e3), B3)

par(mfrow = c(2,2))
for(I in 1:D){
  ylim_ = range(c(c1[,I], c2[,I], c3[,I]))
  plot(ord, c1[,I], type = 'l', ylim = ylim_)
  points(ord, c2[,I], type='l', col = 'red', lty = 2)
  points(ord, c3[,I], type='l', col = 'blue', lty = 3)
}

C1 = cbind(data.frame(c1))
C1$ind = 1:nrow(C1)/nrow(C1)
library(ggtern)
df.lims = data.frame(x1 = c(0,.2,.8),
                     x2 = c(0,.8,.00),
                     x3 = c(1,.00,.2))
ggtern() +
  geom_line(data=C1, aes(x1,x2,x3,col=ind,alpha=ind)) +
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

ggtern() +
  geom_line(data=C1, aes(x1,x2,x4,col=ind,alpha=ind))+
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

# > par('din')
# [1] 8.645833 2.760417
