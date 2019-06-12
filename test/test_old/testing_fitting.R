library(coda.count)
set.seed(1)
n = 10
MU.orig = c(0,0)
SIGMA.orig = diag(3, 2)
X = rlrnm(MU.orig, SIGMA.orig, n, 100)

ft1 = c_lrnm_fit_hermite(X, MU.orig, sigma0 = SIGMA.orig)
ft1[[1]]
ft1[[2]]

library(randtoolbox)
Z = randtoolbox::sobol(n = n, d = 2)
ft2 = c_lrnm_fit_mc(X, MU.orig, sigma0 = SIGMA.orig, Z = Z)
ft2[[1]]
ft2[[2]]

ft2 = fit_lrnm(X, method = 'maximum')
ft2$mu
ft2$sigma

pars.dm = fit_dm(X)
H = ilr_coordinates(t(t(X) + pars.dm[,1]))

fit.lrnm = c_lrnm_fit_hermite(X, colMeans(H), cov(H))
fit.lrnm2 = c_lrnm_fit_hermite(X, colMeans(H), cov(H), max_steps = 0, order = 250)

#####3
NSIM = 100000
Z = sobol(NSIM, dim = length(MU.orig), normal = TRUE, init=TRUE)
fit.lrnm2[1:2]
c_lrnm_fit_mc(X, colMeans(H), cov(H), Z, tol=0.0001)[1:2]
normalmultinomial::nm_fit(X,eps = 0.0001)[1:2]

mm = c_lrnm_fit_maximum(X, colMeans(H), cov(H), tol = 0.0001)[1:2]

