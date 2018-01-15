

##############
library(readxl)
library(dplyr)
library(coda.count)
library(parallel)
library(randtoolbox)

sheets = excel_sheets('~/donato/database_e_carote_conteggi_v2_v5_v6.xls')

fitting = list()
times = list()
for(sheet in sheets){
  data = read_excel('~/donato/database_e_carote_conteggi_v2_v5_v6.xls', sheet = sheet, col_names = FALSE)
  vnames = LETTERS[1:ncol(data)]
  nc = c(1, ncol(data), which(colSums(data) == 0))
  X = data %>%
    dplyr::select(-nc) %>%
    as.matrix() %>%
    unname()
  vnames = vnames[-nc]
  H.dm = ilr_coordinates(t(t(X) + fit_dm(X)[,1]))

  MU = colMeans(H.dm)
  SIGMA = cov(H.dm)

  set.seed(1)
  NSIM = 100000
  Z = sobol(NSIM, dim = ncol(H.dm), normal = TRUE, init=TRUE)

  t.start = proc.time()
  fitting[[sheet]] = c_lrnm_fit_mc(X, MU, SIGMA, Z, tol = 0.01, em_max_steps = 50)
  colnames(fitting[[sheet]][[3]]) = vnames
  times[[sheet]] = proc.time() - t.start
}

save(sheets, fitting, times, file = 'test/fitting2.RData')

#
#   mu = MU
#   sigma = SIGMA
#
#   M1i = array(0, dim = c(ncol(Z),nrow(X)))
#   M2i = array(0, dim = c(ncol(Z),ncol(Z),nrow(X)))
#   M2 = array(0, dim = c(ncol(Z),ncol(Z)))
#   Hs = matrix(0, ncol=ncol(Z), nrow=nrow(Z))
#
#   fit_init = c_lrnm_fit_mc_init(X, mu, sigma, Z)
#
#   M1i = array(fit_init[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_init[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step1 = c_lrnm_fit_mc_step(X, fit_init[[1]][,1], fit_init[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step1[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step1[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step2 = c_lrnm_fit_mc_step(X, fit_step1[[1]][,1], fit_step1[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step2[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step2[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step3 = c_lrnm_fit_mc_step(X, fit_step2[[1]][,1], fit_step2[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step3[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step3[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step4 = c_lrnm_fit_mc_step(X, fit_step3[[1]][,1], fit_step3[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step4[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step4[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step5 = c_lrnm_fit_mc_step(X, fit_step4[[1]][,1], fit_step4[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step5[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step5[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step6 = c_lrnm_fit_mc_step(X, fit_step5[[1]][,1], fit_step5[[2]], Z, M1i, M2i)
#
#   M1i = array(fit_step5[[3]], dim = c(ncol(Z),nrow(X)))
#   M2i = array(fit_step5[[4]], dim = c(ncol(Z),ncol(Z),nrow(X)))
#   fit_step6 = c_lrnm_fit_mc_step(X, fit_step5[[1]][,1], fit_step5[[2]], Z, M1i, M2i)
#
#   #########
#   #######
#   i = 159
#   x = X[i,]
#   sampling_mu = fit_init[[3]][,i]
#   sampling_sigma = fit_init[[4]][,,i]
#
#   lik_std = expected_mc_03_init(x, MU, SIGMA, Z, sampling_mu, sampling_sigma, Hs)
#
#
#   for(i in 1:nrow(X)){
#     x = X[i,]
#
#     mu_sampling = mvf_maximum(x, MU, SIGMA, ilr_basis(ncol(X)), 0.0001, 50, 0.66)
#
#     lik_std = expected_mc_01_init(x, MU, SIGMA, Z, mu_sampling, Hs)
#     M1i[,i]  = expected_mc_mean(x, Hs, lik_std)
#     M2 = M2 + expected_mc_var(x, mu, Hs, lik_std)
#     M2i[,,i] = expected_mc_var(x, M1i[,i], Hs, lik_std)
#   }
#
#
#   t.diff = proc.time() - t.start
#
# #}
