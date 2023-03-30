#include <RcppArmadillo.h>

#ifndef DM_H
#define DM_H

Rcpp::List c_dm_fit(arma::mat& X, double eps, int maxiter);

#endif
