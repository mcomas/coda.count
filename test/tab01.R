SEED = 1
library(coda.dist)
library(randtoolbox)

x = c(1,0)
MU = 0
SIGMA = matrix(1)

NSIM = 1000

f0 = function() expected_hermite(x = x, mu_ilr = MU, sigma_ilr = SIGMA, order = 2000)
r0 = f0()



f1 = function(){
  set.seed(SEED)
  fit = function(){
    Z = rnorm(NSIM*length(MU))
    Z = matrix(Z, ncol=length(MU))
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}

f2 = function(){
  set.seed(SEED)
  fit = function(){
    Z0 = rnorm( (NSIM/2)*length(MU))
    Z = matrix(0, nrow=NSIM, ncol=1)
    Z[seq(1,NSIM,2),] = Z0
    Z[seq(2,NSIM,2),] = -Z0
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}

f3 = function(){
  set.seed(SEED)
  fit = function(){
    Z = halton(NSIM, dim = length(MU), normal = TRUE, init = FALSE)
    Z = matrix(Z, ncol=length(MU))
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}
f4 = function(){
  set.seed(SEED)
  fit = function(){
    Z0 = halton(NSIM/2, dim = length(MU), normal = TRUE, init=FALSE)
    Z = matrix(0, nrow=NSIM, ncol=1)
    Z[seq(1,NSIM,2),] = Z0
    Z[seq(2,NSIM,2),] = -Z0
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}
f3b = function(){
  set.seed(SEED)
  fit = function(){
    Z = sobol(NSIM, dim = length(MU), normal = TRUE, scrambling = 2, seed = as.integer(runif(1) * 10000))
    Z = matrix(Z, ncol=length(MU))
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}
f4b = function(){
  set.seed(SEED)
  fit = function(){
    Z0 = sobol(NSIM/2, dim = length(MU), normal = TRUE, scrambling = 2, seed = as.integer(runif(1) * 10000))
    Z = matrix(0, nrow=NSIM, ncol=1)
    Z[seq(1,NSIM,2),] = Z0
    Z[seq(2,NSIM,2),] = -Z0
    E = expected_montecarlo_01(x = x, mu_ilr = MU, sigma_ilr = SIGMA,
                               Z = Z, mu_exp = r0[[2]])
    c(E[[2]][1],E[[3]][1])
  }

  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))
}

f5 = function(){
  set.seed(SEED)
  fit = function(){
    E = expected_metropolis(x = x, mu_ilr = MU, sigma_ilr = SIGMA, mu_exp = r0[[2]], nsim = NSIM)
    c(E[[2]][1],E[[3]][1])
  }
  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))

}

library(rstan)
load('test/nm_stan_model.RData')
# MCMC using hamiltonian algorithm
f6 = function(){
  set.seed(SEED)
  fit = function(){
    E = expected_hamiltonian(x = x, mu_ilr = MU, sigma_ilr = SIGMA, nsim = NSIM)
    c(E[[2]][1],E[[3]][1])
  }
  res = t(replicate(500, fit()))
  list('m' = apply(res, 2, mean),
       'sd' = apply(res, 2, sd))

}

r1 = f1()
r2 = f2()
r3 = f3()
r4 = f4()
r5 = f5()
r6 = f6()

