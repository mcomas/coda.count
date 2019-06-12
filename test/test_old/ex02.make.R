library(dcount)
library(rstan)
setwd('test')

load('stan_model.RData')
KS = 3
SEEDS = 1:500

COMB = expand.grid(KS, SEEDS)
for(i in 1:NROW(COMB)){
  row = as.numeric(COMB[i,])
  PATTERN = sprintf('K_%04d-SEED_%04d', row[1], row[2])
  print(PATTERN)
  source('ex02.R')
  save(r0,r1,r2,r3,r4,r5,r6, file = sprintf('ex02_R/data-%s.RData', PATTERN))
  list(r0,r1,r2,r3,r4,r5,r6)
  rm(PATTERN)
}

