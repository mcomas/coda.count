# include <RcppArmadillo.h>

#ifndef CODA_BASE_H
#define CODA_BASE_H

arma::mat inv_B_coordinates(arma::mat ilrX, arma::mat Binv);
arma::mat ilr_basis_default(unsigned int dim);

#endif
