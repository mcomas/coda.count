library(coda.count)
library(coda.base)
load('data/parliament2015.RData')

Y_all = as.matrix(parliament2015[,c('jxsi', 'psc', 'pp', 'catsp', 'cs', 'cup')])

Y = with(parliament2015, cbind('cat' = jxsi+cup, 'esp' = pp+psc+cs, catsp))
X = model.matrix(~I(log(pop))+I(100*birth.cat/pop), data = parliament2015)

alpha = c_dm_fit_alpha(Y)
P = t(t(Y) + alpha[,1])
H = coordinates(P)

# m = lm(H~0+X)
# b = composition(coef(m))
# b/b[,1]

D = ncol(Y)
B = ilr_basis(D)

# c_fit_lrnm_hermite(Y, B, order = 5)
# c_fit_lrnm_hermite_precision(Y, B, order = 5)
# c_fit_lm_lrnm_hermite(Y, B, X[,1:3, drop=FALSE], order = 5)
fit = c_fit_lm_lrnm_hermite_centered(Y, B, X[,1:3, drop=FALSE], order = 5)
fitted.values = X %*% fit[[1]]
P.fitted.values = composition(fitted.values, B)
colnames(P.fitted.values) = colnames(Y)
library(ggtern)

plot(X[,2], fitted.values[,1])

fit[[1]]
X[91,]
scatter3D(x = X[,2], y = X[,3], z = fitted.values[,1], theta = 50, phi = 0)

ggtern(data = as.data.frame(Y)) +
  geom_point(aes(x = cat, y = esp, z = catsp), size = 2) +
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

ggtern(data = as.data.frame(P.fitted.values)) +
  geom_point(aes(x = cat, y = esp, z = catsp), size = 2) +
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

#
# head(Y)
#
# m = lm(Y~0+X)
#
# coef(m)
# solve(t(X) %*% X) %*% t(X) %*% Y
#
# vcov(m)[1:3,1:3]
# res = (Y-X %*% coef(m))
# S_res = t(res) %*% res / (nrow(X)-3)
# kronecker(solve(t(X) %*% X),  S_res)[c(1,4,7), c(1,4,7)]
# # arma::kron( A, B )
#
# res = (Y-X %*% coef(m))
# head(t(res) %*% res)
# head(res)
# head(m$residuals)
# cov(m$residuals) *(nrow(X)-1)/(nrow(X)-3)
# (t(res) %*% res) / (nrow(X)-3)
#
#
# inc_mu = function(n, x, mu){
#   mu + (x - mu)/n
# }
# inc_sigma = function(n, x, mu, sigma){
#   sigma + (x-mu)^2/n - sigma/(n-1)
# }
# inc_sigma_n = function(n, x, mu, mu_prev, sigma_n){
#   sigma_n + ((x-mu)*(x-mu_prev)-sigma_n)/n
# }
# (mu = inc_mu(1, 20, 0))
#
# (sigma = inc_sigma(2, 15, mu, 0))
# (mu = inc_mu(2, 15, mu))
# (sigma_n = inc_sigma_n(2, 15, mu, 0))
#
# (sigma = inc_sigma(3, 22, mu, sigma))
# (mu = inc_mu(3, 22, mu))
#
