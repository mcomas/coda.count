SEED = 1000
library(coda.count)
library(randtoolbox)

x = c(1,0)
MU = 0
SIGMA = matrix(1)

NSIM = 1000

f0 = function() expected_hermite(x = x, mu_ilr = MU, sigma_ilr = SIGMA, order = 2000)
r0 = f0()

step_by_step = function(Z){
  res = sapply(1:nrow(Z), function(i){
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = matrix(Z[1:i,],ncol=ncol(Z)), mu_exp = r0[[2]])
    c(E[[2]][1], E[[3]][1])
  })
  t(res)
}


f1 = function(){
  set.seed(SEED)
  Z = rnorm(NSIM*length(MU))
  Z = matrix(Z, ncol=length(MU))
  step_by_step(Z)
}
f2 = function(){
  set.seed(SEED)
  Z0 = rnorm( (NSIM/2)*length(MU))
  Z = matrix(0, nrow=NSIM, ncol=1)
  Z[seq(1,NSIM,2),] = Z0
  Z[seq(2,NSIM,2),] = -Z0
  step_by_step(Z)
}
f3 = function(){
  set.seed(SEED)
  Z = halton(NSIM, dim = length(MU), normal = TRUE)
  Z = matrix(Z, ncol=length(MU))
  step_by_step(Z)
}
f4 = function(){
  set.seed(SEED)
  Z0 = halton(NSIM/2, dim = length(MU), normal = TRUE)
  Z = matrix(0, nrow=NSIM, ncol=1)
  Z[seq(1,NSIM,2),] = Z0
  Z[seq(2,NSIM,2),] = -Z0
  step_by_step(Z)
}
f5 = function(){
  res = sapply(1:NSIM, function(i){
    set.seed(SEED)
    E = expected_metropolis(x = x, mu_ilr = MU, sigma_ilr = SIGMA, mu_exp = r0[[2]], nsim = i)
    c(E[[2]][1], E[[3]][1])
  })
  t(res)
}

library(rstan)
load('test/nm_stan_model.RData')
# MCMC using hamiltonian algorithm
f6 = function(){
  res = sapply(1:NSIM, function(i){
    set.seed(SEED)
    E = expected_hamiltonian(x = x, mu_ilr = MU, sigma_ilr = SIGMA, nsim = i)
    c(E[[2]][1], E[[3]][1])
  })
  t(res)
}

r1 = f1()
r2 = f2()
r3 = f3()
r4 = f4()
r5 = f5()
r6 = f6()

XLIM = 200:1000
YLIM = range(r1[XLIM,1], r2[XLIM,1], r3[XLIM,1],
             r4[XLIM,1], r5[XLIM,1], r6[XLIM,1])
plot(r1[,1], type='l', xlim = range(XLIM), ylim = YLIM,
     xlab = 'Iteration', ylab = 'Estimation')
segments(x0 = 1, x1 = NSIM, y0 = r0[[2]], y1 = r0[[2]], col='red')
points(r2[,1], col='blue', type='l')
points(r3[,1], col='green', type='l')
points(r4[,1], col='orange', type='l')
points(r5[,1], col='brown', type='l')
points(r6[,1], col='purple', type='l')
legend('bottomright',
       legend = c('MC', 'MC (Antithetic variates)','QMC', 'QMC (Antithetic variates)',
                  'MCMC (Metropolis)', 'MCMC (Hamiltonian)'),
       col= c('black', 'blue', 'green', 'orange', 'brown', 'purple'),
       bty = 'n', lty=1, cex = 0.8)

# WIDTH = 6.5, HEIGHT = 5.5
XLIM = 200:1000
YLIM = range(r1[XLIM,2], r2[XLIM,2], r3[XLIM,2],
             r4[XLIM,2], r5[XLIM,2], r6[XLIM,2])
plot(r1[,2], type='l', xlim = range(XLIM), ylim = YLIM,
     xlab = 'Iteration', ylab = 'Estimation')
segments(x0 = 1, x1 = NSIM, y0 = r0[[3]], y1 = r0[[3]], col='red')
points(r2[,2], col='blue', type='l')
points(r3[,2], col='green', type='l')
points(r4[,2], col='orange', type='l')
points(r5[,2], col='brown', type='l')
points(r6[,2], col='purple', type='l')
legend('bottomright',
       legend = c('MC', 'MC (Antithetic variates)','QMC', 'QMC (Antithetic variates)',
                  'MCMC (Metropolis)', 'MCMC (Hamiltonian)'),
       col= c('black', 'blue', 'green', 'orange', 'brown', 'purple'),
       bty = 'n', lty=1, cex = 0.8)

