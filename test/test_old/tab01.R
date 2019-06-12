SEED = 1
library(coda.dist)
library(randtoolbox)
library(dplyr)
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

add_time = function(f){
  t0 = proc.time()
  res = f()
  attr(res, 'time') = proc.time() - t0
  res
}
r1 = add_time(f1)
r2 = add_time(f2)
r3 = add_time(f3)
r4 = add_time(f4)
r5 = add_time(f5)
r6 = add_time(f6)

r = list(r1, r2, r3, r4, r5, r6)
M = sapply(r, nth, 1)
V = sapply(r, nth, 2)
TIME = sapply(r, attr, 'time')[3,]
TIME = TIME / TIME[1]
tab.out = data.frame(
  'Method' = c('MC', 'MC (Antithetic variates)','QMC', 'QMC (Antithetic variates)',
               'MCMC (Metropolis)', 'MCMC (Hamiltonian)'),
  'First moment' = sprintf("%8.7f (%6.5f)", M[1,], V[1,]),
  'Second momenth' = sprintf("%8.7f (%6.5f)", M[2,], V[2,]),
  'Time' = sprintf("x%.2f", TIME))
tab.out

# sink(file='tex/convergence01-quasi.tex')
cat(Hmisc::latexTabular(tab.out,
                        headings = sprintf("\\textbf{%s}", names(tab.out)),
                        hline = 1, align = 'l | r r | l', helvetica = F, translate = F))
# sink()
