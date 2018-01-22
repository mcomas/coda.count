library(coda.count)

MU.orig = c(0,0)
SIGMA.orig = diag(2)
X = rlrnm(MU.orig, SIGMA.orig, 100, 10)

fit_lrnm(X)

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

