
#ifndef GAUSSIAN_APPROX_H
#define GAUSSIAN_APPROX_H

arma::mat c_dirichlet_alr_approx(arma::vec alpha);

arma::mat c_logistic_ilr_approximation(int d);

arma::mat c_logistic_alr_approximation(int d);

arma::mat c_gaussian_product(arma::mat pars1, arma::mat pars2);

arma::mat c_gaussian_division(arma::mat pars1, arma::mat pars_res);

#endif
