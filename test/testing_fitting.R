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
NSIM = 1000
Z0 = sobol(NSIM/2, dim = length(MU.orig), normal = TRUE, init=TRUE)
Z = matrix(0, nrow=NSIM, ncol=length(MU.orig))
Z[seq(1,NSIM,2),] = Z0
Z[seq(2,NSIM,2),] = -Z0

c_lrnm_fit_mc(X, colMeans(H), cov(H), Z)[1:2]
normalmultinomial::nm_fit(X)[1:2]

mm = c_lrnm_fit_maximum(X, colMeans(H), cov(H))[1:2]

