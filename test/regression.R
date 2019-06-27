library(coda.count)
load('data/parliament2015.RData')

Y_all = as.matrix(parliament2015[,c('jxsi', 'psc', 'pp', 'catsp', 'cs', 'cup')])

Y = with(parliament2015, cbind('cat' = jxsi+cup, 'esp' = pp+psc+cs, catsp))
X = model.matrix(~pop+birth.cat, data = parliament2015)

D = ncol(Y)
B = ilr_basis(D)

c_fit_lrnm_hermite(Y, B, order = 5)
c_fit_lrnm_hermite_precision(Y, B, order = 5)
c_fit_lm_lrnm_hermite(Y, B, X[,1, drop=FALSE], order = 5)
