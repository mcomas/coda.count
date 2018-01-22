// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ilr_basis
arma::mat ilr_basis(unsigned int dim);
RcppExport SEXP _coda_count_ilr_basis(SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(ilr_basis(dim));
    return rcpp_result_gen;
END_RCPP
}
// ilr_coordinates
arma::mat ilr_coordinates(arma::mat X);
RcppExport SEXP _coda_count_ilr_coordinates(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ilr_coordinates(X));
    return rcpp_result_gen;
END_RCPP
}
// inv_ilr_coordinates
arma::mat inv_ilr_coordinates(arma::mat ilrX);
RcppExport SEXP _coda_count_inv_ilr_coordinates(SEXP ilrXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ilrX(ilrXSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_ilr_coordinates(ilrX));
    return rcpp_result_gen;
END_RCPP
}
// ilr_to_alr
arma::mat ilr_to_alr(arma::mat ILR);
RcppExport SEXP _coda_count_ilr_to_alr(SEXP ILRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ILR(ILRSEXP);
    rcpp_result_gen = Rcpp::wrap(ilr_to_alr(ILR));
    return rcpp_result_gen;
END_RCPP
}
// alr_to_ilr
arma::mat alr_to_ilr(arma::mat ALR);
RcppExport SEXP _coda_count_alr_to_ilr(SEXP ALRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ALR(ALRSEXP);
    rcpp_result_gen = Rcpp::wrap(alr_to_ilr(ALR));
    return rcpp_result_gen;
END_RCPP
}
// c_ddm
double c_ddm(arma::vec x, arma::vec alpha);
RcppExport SEXP _coda_count_c_ddm(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_ddm(x, alpha));
    return rcpp_result_gen;
END_RCPP
}
// c_dm_fit
Rcpp::List c_dm_fit(arma::mat X, double eps, int maxiter);
RcppExport SEXP _coda_count_c_dm_fit(SEXP XSEXP, SEXP epsSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dm_fit(X, eps, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _coda_count_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// validate_dnm
double validate_dnm(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);
RcppExport SEXP _coda_count_validate_dnm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(validate_dnm(x, mu, sigma, order));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_01
Rcpp::List expected_montecarlo_01(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_exp);
RcppExport SEXP _coda_count_expected_montecarlo_01(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_exp(mu_expSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_01(x, mu_ilr, sigma_ilr, Z, mu_exp));
    return rcpp_result_gen;
END_RCPP
}
// expected_mc_01_init
arma::vec expected_mc_01_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& sigma_ilr, arma::mat& Z, arma::vec& sampling_mu, arma::mat& Hs);
RcppExport SEXP _coda_count_expected_mc_01_init(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP sampling_muSEXP, SEXP HsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sampling_mu(sampling_muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hs(HsSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_mc_01_init(x, mu_ilr, sigma_ilr, Z, sampling_mu, Hs));
    return rcpp_result_gen;
END_RCPP
}
// expected_mc_03_init
arma::vec expected_mc_03_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& inv_sigma_ilr, arma::mat& Z, arma::vec& mu_sampling, arma::mat& sigma_sampling, arma::mat& Hs);
RcppExport SEXP _coda_count_expected_mc_03_init(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP inv_sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_samplingSEXP, SEXP sigma_samplingSEXP, SEXP HsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma_ilr(inv_sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu_sampling(mu_samplingSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma_sampling(sigma_samplingSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hs(HsSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_mc_03_init(x, mu_ilr, inv_sigma_ilr, Z, mu_sampling, sigma_sampling, Hs));
    return rcpp_result_gen;
END_RCPP
}
// expected_mc_mean
arma::vec expected_mc_mean(arma::vec& x, arma::mat& Hs, arma::vec& lik_st);
RcppExport SEXP _coda_count_expected_mc_mean(SEXP xSEXP, SEXP HsSEXP, SEXP lik_stSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hs(HsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lik_st(lik_stSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_mc_mean(x, Hs, lik_st));
    return rcpp_result_gen;
END_RCPP
}
// expected_mc_var
arma::mat expected_mc_var(arma::vec& x, arma::vec& mu, arma::mat& Hs, arma::vec& lik_st);
RcppExport SEXP _coda_count_expected_mc_var(SEXP xSEXP, SEXP muSEXP, SEXP HsSEXP, SEXP lik_stSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hs(HsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lik_st(lik_stSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_mc_var(x, mu, Hs, lik_st));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_01_mean
arma::vec expected_montecarlo_01_mean(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_exp);
RcppExport SEXP _coda_count_expected_montecarlo_01_mean(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_exp(mu_expSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_01_mean(x, mu_ilr, sigma_ilr, Z, mu_exp));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_01_centered
Rcpp::List expected_montecarlo_01_centered(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_exp);
RcppExport SEXP _coda_count_expected_montecarlo_01_centered(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_exp(mu_expSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_01_centered(x, mu_ilr, sigma_ilr, Z, mu_exp));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_02
arma::mat expected_montecarlo_02(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_exp, double var_exp);
RcppExport SEXP _coda_count_expected_montecarlo_02(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_expSEXP, SEXP var_expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_exp(mu_expSEXP);
    Rcpp::traits::input_parameter< double >::type var_exp(var_expSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_02(x, mu_ilr, sigma_ilr, Z, mu_exp, var_exp));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_03
Rcpp::List expected_montecarlo_03(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_sampling, arma::mat sigma_sampling);
RcppExport SEXP _coda_count_expected_montecarlo_03(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_samplingSEXP, SEXP sigma_samplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_sampling(mu_samplingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_sampling(sigma_samplingSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_03(x, mu_ilr, sigma_ilr, Z, mu_sampling, sigma_sampling));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo_04
Rcpp::List expected_montecarlo_04(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec m1, arma::mat m2);
RcppExport SEXP _coda_count_expected_montecarlo_04(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo_04(x, mu_ilr, sigma_ilr, Z, m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// expected_hermite
Rcpp::List expected_hermite(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, int order);
RcppExport SEXP _coda_count_expected_hermite(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_hermite(x, mu_ilr, sigma_ilr, order));
    return rcpp_result_gen;
END_RCPP
}
// expected_montecarlo
Rcpp::List expected_montecarlo(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat& Z, arma::vec mu_sampling, arma::mat sigma_sampling, arma::mat& Hz);
RcppExport SEXP _coda_count_expected_montecarlo(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP ZSEXP, SEXP mu_samplingSEXP, SEXP sigma_samplingSEXP, SEXP HzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_sampling(mu_samplingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_sampling(sigma_samplingSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hz(HzSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_montecarlo(x, mu_ilr, sigma_ilr, Z, mu_sampling, sigma_sampling, Hz));
    return rcpp_result_gen;
END_RCPP
}
// expected_metropolis
Rcpp::List expected_metropolis(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::vec mu_exp, int nsim, int ignored_steps);
RcppExport SEXP _coda_count_expected_metropolis(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP mu_expSEXP, SEXP nsimSEXP, SEXP ignored_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_exp(mu_expSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type ignored_steps(ignored_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_metropolis(x, mu_ilr, sigma_ilr, mu_exp, nsim, ignored_steps));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_sample
arma::mat metropolis_sample(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::vec x0, int nsim, int ignored_steps);
RcppExport SEXP _coda_count_metropolis_sample(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP x0SEXP, SEXP nsimSEXP, SEXP ignored_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type ignored_steps(ignored_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_sample(x, mu_ilr, sigma_ilr, x0, nsim, ignored_steps));
    return rcpp_result_gen;
END_RCPP
}
// mvf_deriv
double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP _coda_count_mvf_deriv(SEXP ISEXP, SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mvf_deriv(I, a, mu, inv_sigma, x));
    return rcpp_result_gen;
END_RCPP
}
// mvf_deriv2
double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP _coda_count_mvf_deriv2(SEXP ISEXP, SEXP JSEXP, SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mvf_deriv2(I, J, a, mu, inv_sigma, x));
    return rcpp_result_gen;
END_RCPP
}
// alr_basis
arma::mat alr_basis(int K);
RcppExport SEXP _coda_count_alr_basis(SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(alr_basis(K));
    return rcpp_result_gen;
END_RCPP
}
// mvf_maximum
arma::vec mvf_maximum(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat B, double eps, int max_iter, double prop);
RcppExport SEXP _coda_count_mvf_maximum(SEXP xSEXP, SEXP mu_ilrSEXP, SEXP sigma_ilrSEXP, SEXP BSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_ilr(mu_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_ilr(sigma_ilrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    rcpp_result_gen = Rcpp::wrap(mvf_maximum(x, mu_ilr, sigma_ilr, B, eps, max_iter, prop));
    return rcpp_result_gen;
END_RCPP
}
// mvf_maximum_alr
arma::vec mvf_maximum_alr(arma::vec x, arma::vec mu_alr, arma::mat& inv_sigma_alr, arma::vec a, double eps, int max_iter);
RcppExport SEXP _coda_count_mvf_maximum_alr(SEXP xSEXP, SEXP mu_alrSEXP, SEXP inv_sigma_alrSEXP, SEXP aSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_alr(mu_alrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma_alr(inv_sigma_alrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(mvf_maximum_alr(x, mu_alr, inv_sigma_alr, a, eps, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// ldnormal
arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP _coda_count_ldnormal(SEXP HSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ldnormal(H, mu, inv_sigma));
    return rcpp_result_gen;
END_RCPP
}
// lpmultinomial_const
double lpmultinomial_const(arma::vec x);
RcppExport SEXP _coda_count_lpmultinomial_const(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(lpmultinomial_const(x));
    return rcpp_result_gen;
END_RCPP
}
// c_dlrnm_hermite
double c_dlrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, int order, int step_by, double eps, int max_steps);
RcppExport SEXP _coda_count_c_dlrnm_hermite(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type step_by(step_bySEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dlrnm_hermite(x, mu, sigma, order, step_by, eps, max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_hermite
Rcpp::List c_lrnm_fit_hermite(arma::mat X, arma::vec mu0, arma::mat sigma0, int order, int step_by, double eps, int max_steps, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_hermite(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type step_by(step_bySEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_hermite(X, mu0, sigma0, order, step_by, eps, max_steps, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_maximum
Rcpp::List c_lrnm_fit_maximum(arma::mat X, arma::vec mu0, arma::mat sigma0, double tol, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_maximum(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP tolSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_maximum(X, mu0, sigma0, tol, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_maximum_alr
Rcpp::List c_lrnm_fit_maximum_alr(arma::mat X, arma::vec mu0, arma::mat sigma0, double tol, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_maximum_alr(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP tolSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_maximum_alr(X, mu0, sigma0, tol, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_mc_init
Rcpp::List c_lrnm_fit_mc_init(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z, double tol, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_mc_init(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP ZSEXP, SEXP tolSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_mc_init(X, mu0, sigma0, Z, tol, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_mc_step
Rcpp::List c_lrnm_fit_mc_step(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z, arma::mat M1i, arma::cube M2i, double tol, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_mc_step(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP ZSEXP, SEXP M1iSEXP, SEXP M2iSEXP, SEXP tolSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M1i(M1iSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type M2i(M2iSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_mc_step(X, mu0, sigma0, Z, M1i, M2i, tol, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_mc
Rcpp::List c_lrnm_fit_mc(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z, double tol, int em_max_steps);
RcppExport SEXP _coda_count_c_lrnm_fit_mc(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP ZSEXP, SEXP tolSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_mc(X, mu0, sigma0, Z, tol, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_rmultinomial_Rcpp
IntegerMatrix c_rmultinomial_Rcpp(NumericMatrix P, NumericVector vsize);
RcppExport SEXP _coda_count_c_rmultinomial_Rcpp(SEXP PSEXP, SEXP vsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vsize(vsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rmultinomial_Rcpp(P, vsize));
    return rcpp_result_gen;
END_RCPP
}
// c_rmultinomial
arma::mat c_rmultinomial(arma::mat P, arma::vec vsize);
RcppExport SEXP _coda_count_c_rmultinomial(SEXP PSEXP, SEXP vsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vsize(vsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rmultinomial(P, vsize));
    return rcpp_result_gen;
END_RCPP
}
// c_rdirichlet
arma::mat c_rdirichlet(int n, arma::vec alpha);
RcppExport SEXP _coda_count_c_rdirichlet(SEXP nSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rdirichlet(n, alpha));
    return rcpp_result_gen;
END_RCPP
}
// c_rdirichletmultinomial
List c_rdirichletmultinomial(arma::vec alpha, arma::vec size);
RcppExport SEXP _coda_count_c_rdirichletmultinomial(SEXP alphaSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rdirichletmultinomial(alpha, size));
    return rcpp_result_gen;
END_RCPP
}
// c_rnormal
arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _coda_count_c_rnormal(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rnormal(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// c_rnormalSimplex
arma::mat c_rnormalSimplex(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _coda_count_c_rnormalSimplex(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rnormalSimplex(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// c_rnormalmultinomial
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size);
RcppExport SEXP _coda_count_c_rnormalmultinomial(SEXP muSEXP, SEXP sigmaSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rnormalmultinomial(mu, sigma, size));
    return rcpp_result_gen;
END_RCPP
}
// c_simplex_lattice
NumericMatrix c_simplex_lattice(int K, int SIZE);
RcppExport SEXP _coda_count_c_simplex_lattice(SEXP KSEXP, SEXP SIZESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type SIZE(SIZESEXP);
    rcpp_result_gen = Rcpp::wrap(c_simplex_lattice(K, SIZE));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_coda_count_ilr_basis", (DL_FUNC) &_coda_count_ilr_basis, 1},
    {"_coda_count_ilr_coordinates", (DL_FUNC) &_coda_count_ilr_coordinates, 1},
    {"_coda_count_inv_ilr_coordinates", (DL_FUNC) &_coda_count_inv_ilr_coordinates, 1},
    {"_coda_count_ilr_to_alr", (DL_FUNC) &_coda_count_ilr_to_alr, 1},
    {"_coda_count_alr_to_ilr", (DL_FUNC) &_coda_count_alr_to_ilr, 1},
    {"_coda_count_c_ddm", (DL_FUNC) &_coda_count_c_ddm, 2},
    {"_coda_count_c_dm_fit", (DL_FUNC) &_coda_count_c_dm_fit, 3},
    {"_coda_count_timesTwo", (DL_FUNC) &_coda_count_timesTwo, 1},
    {"_coda_count_validate_dnm", (DL_FUNC) &_coda_count_validate_dnm, 4},
    {"_coda_count_expected_montecarlo_01", (DL_FUNC) &_coda_count_expected_montecarlo_01, 5},
    {"_coda_count_expected_mc_01_init", (DL_FUNC) &_coda_count_expected_mc_01_init, 6},
    {"_coda_count_expected_mc_03_init", (DL_FUNC) &_coda_count_expected_mc_03_init, 7},
    {"_coda_count_expected_mc_mean", (DL_FUNC) &_coda_count_expected_mc_mean, 3},
    {"_coda_count_expected_mc_var", (DL_FUNC) &_coda_count_expected_mc_var, 4},
    {"_coda_count_expected_montecarlo_01_mean", (DL_FUNC) &_coda_count_expected_montecarlo_01_mean, 5},
    {"_coda_count_expected_montecarlo_01_centered", (DL_FUNC) &_coda_count_expected_montecarlo_01_centered, 5},
    {"_coda_count_expected_montecarlo_02", (DL_FUNC) &_coda_count_expected_montecarlo_02, 6},
    {"_coda_count_expected_montecarlo_03", (DL_FUNC) &_coda_count_expected_montecarlo_03, 6},
    {"_coda_count_expected_montecarlo_04", (DL_FUNC) &_coda_count_expected_montecarlo_04, 6},
    {"_coda_count_expected_hermite", (DL_FUNC) &_coda_count_expected_hermite, 4},
    {"_coda_count_expected_montecarlo", (DL_FUNC) &_coda_count_expected_montecarlo, 7},
    {"_coda_count_expected_metropolis", (DL_FUNC) &_coda_count_expected_metropolis, 6},
    {"_coda_count_metropolis_sample", (DL_FUNC) &_coda_count_metropolis_sample, 6},
    {"_coda_count_mvf_deriv", (DL_FUNC) &_coda_count_mvf_deriv, 5},
    {"_coda_count_mvf_deriv2", (DL_FUNC) &_coda_count_mvf_deriv2, 6},
    {"_coda_count_alr_basis", (DL_FUNC) &_coda_count_alr_basis, 1},
    {"_coda_count_mvf_maximum", (DL_FUNC) &_coda_count_mvf_maximum, 7},
    {"_coda_count_mvf_maximum_alr", (DL_FUNC) &_coda_count_mvf_maximum_alr, 6},
    {"_coda_count_ldnormal", (DL_FUNC) &_coda_count_ldnormal, 3},
    {"_coda_count_lpmultinomial_const", (DL_FUNC) &_coda_count_lpmultinomial_const, 1},
    {"_coda_count_c_dlrnm_hermite", (DL_FUNC) &_coda_count_c_dlrnm_hermite, 7},
    {"_coda_count_c_lrnm_fit_hermite", (DL_FUNC) &_coda_count_c_lrnm_fit_hermite, 8},
    {"_coda_count_c_lrnm_fit_maximum", (DL_FUNC) &_coda_count_c_lrnm_fit_maximum, 5},
    {"_coda_count_c_lrnm_fit_maximum_alr", (DL_FUNC) &_coda_count_c_lrnm_fit_maximum_alr, 5},
    {"_coda_count_c_lrnm_fit_mc_init", (DL_FUNC) &_coda_count_c_lrnm_fit_mc_init, 6},
    {"_coda_count_c_lrnm_fit_mc_step", (DL_FUNC) &_coda_count_c_lrnm_fit_mc_step, 8},
    {"_coda_count_c_lrnm_fit_mc", (DL_FUNC) &_coda_count_c_lrnm_fit_mc, 6},
    {"_coda_count_c_rmultinomial_Rcpp", (DL_FUNC) &_coda_count_c_rmultinomial_Rcpp, 2},
    {"_coda_count_c_rmultinomial", (DL_FUNC) &_coda_count_c_rmultinomial, 2},
    {"_coda_count_c_rdirichlet", (DL_FUNC) &_coda_count_c_rdirichlet, 2},
    {"_coda_count_c_rdirichletmultinomial", (DL_FUNC) &_coda_count_c_rdirichletmultinomial, 2},
    {"_coda_count_c_rnormal", (DL_FUNC) &_coda_count_c_rnormal, 3},
    {"_coda_count_c_rnormalSimplex", (DL_FUNC) &_coda_count_c_rnormalSimplex, 3},
    {"_coda_count_c_rnormalmultinomial", (DL_FUNC) &_coda_count_c_rnormalmultinomial, 3},
    {"_coda_count_c_simplex_lattice", (DL_FUNC) &_coda_count_c_simplex_lattice, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_coda_count(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
