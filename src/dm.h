#include <RcppArmadillo.h>

#ifndef DM_H
#define DM_H

arma::vec c_dm_fit_alpha(arma::mat X, double eps = 0.0001, int maxiter = 5000);

#endif
