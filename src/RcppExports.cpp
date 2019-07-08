// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// H
arma::mat H(int d);
RcppExport SEXP _coda_count_H(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(H(d));
    return rcpp_result_gen;
END_RCPP
}
// F
arma::mat F(int d);
RcppExport SEXP _coda_count_F(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(F(d));
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
// c_dirichlet_alr_approx
arma::mat c_dirichlet_alr_approx(arma::vec alpha);
RcppExport SEXP _coda_count_c_dirichlet_alr_approx(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dirichlet_alr_approx(alpha));
    return rcpp_result_gen;
END_RCPP
}
// c_logistic_ilr_approximation
arma::mat c_logistic_ilr_approximation(int d);
RcppExport SEXP _coda_count_c_logistic_ilr_approximation(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(c_logistic_ilr_approximation(d));
    return rcpp_result_gen;
END_RCPP
}
// c_logistic_alr_approximation
arma::mat c_logistic_alr_approximation(int d);
RcppExport SEXP _coda_count_c_logistic_alr_approximation(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(c_logistic_alr_approximation(d));
    return rcpp_result_gen;
END_RCPP
}
// c_gaussian_product
arma::mat c_gaussian_product(arma::mat pars1, arma::mat pars2);
RcppExport SEXP _coda_count_c_gaussian_product(SEXP pars1SEXP, SEXP pars2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pars1(pars1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pars2(pars2SEXP);
    rcpp_result_gen = Rcpp::wrap(c_gaussian_product(pars1, pars2));
    return rcpp_result_gen;
END_RCPP
}
// c_gaussian_division
arma::mat c_gaussian_division(arma::mat pars1, arma::mat pars_res);
RcppExport SEXP _coda_count_c_gaussian_division(SEXP pars1SEXP, SEXP pars_resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pars1(pars1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pars_res(pars_resSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gaussian_division(pars1, pars_res));
    return rcpp_result_gen;
END_RCPP
}
// c_posterior_approximation
arma::cube c_posterior_approximation(arma::mat X, arma::vec mu, arma::mat& sigma, arma::mat& B);
RcppExport SEXP _coda_count_c_posterior_approximation(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(c_posterior_approximation(X, mu, sigma, B));
    return rcpp_result_gen;
END_RCPP
}
// c_fit_lrnm_gaussian_approx
Rcpp::List c_fit_lrnm_gaussian_approx(arma::mat X, arma::mat B, double eps, int max_iter, int em_max_steps);
RcppExport SEXP _coda_count_c_fit_lrnm_gaussian_approx(SEXP XSEXP, SEXP BSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP em_max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_steps(em_max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_fit_lrnm_gaussian_approx(X, B, eps, max_iter, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_d_lrnm_hermite
double c_d_lrnm_hermite(arma::vec x, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_d_lrnm_hermite(SEXP xSEXP, SEXP mu_priorSEXP, SEXP sigma_priorSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior(mu_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_hermite(x, mu_prior, sigma_prior, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_moments_lrnm_hermite_precision_lm
arma::mat c_moments_lrnm_hermite_precision_lm(arma::vec x, arma::vec mu, arma::mat sigma, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_moments_lrnm_hermite_precision_lm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_priorSEXP, SEXP sigma_priorSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior(mu_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_moments_lrnm_hermite_precision_lm(x, mu, sigma, mu_prior, sigma_prior, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_obtain_moments_lrnm_hermite
Rcpp::List c_obtain_moments_lrnm_hermite(arma::mat Y, arma::vec mu, arma::mat sigma, arma::mat B, int order);
RcppExport SEXP _coda_count_c_obtain_moments_lrnm_hermite(SEXP YSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_obtain_moments_lrnm_hermite(Y, mu, sigma, B, order));
    return rcpp_result_gen;
END_RCPP
}
// c_fit_lm_lrnm_hermite_centered
Rcpp::List c_fit_lm_lrnm_hermite_centered(arma::mat Y, arma::mat B, arma::mat X, int order, double eps, int max_iter);
RcppExport SEXP _coda_count_c_fit_lm_lrnm_hermite_centered(SEXP YSEXP, SEXP BSEXP, SEXP XSEXP, SEXP orderSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_fit_lm_lrnm_hermite_centered(Y, B, X, order, eps, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// c_d_lrnm_montecarlo
double c_d_lrnm_montecarlo(arma::vec x, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv, arma::mat& Z);
RcppExport SEXP _coda_count_c_d_lrnm_montecarlo(SEXP xSEXP, SEXP mu_priorSEXP, SEXP sigma_priorSEXP, SEXP BinvSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior(mu_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_montecarlo(x, mu_prior, sigma_prior, Binv, Z));
    return rcpp_result_gen;
END_RCPP
}
// c_moments_lrnm_montecarlo_precision_lm
arma::mat c_moments_lrnm_montecarlo_precision_lm(arma::vec x, arma::vec mu, arma::mat sigma, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv, arma::mat& Z);
RcppExport SEXP _coda_count_c_moments_lrnm_montecarlo_precision_lm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_priorSEXP, SEXP sigma_priorSEXP, SEXP BinvSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior(mu_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(c_moments_lrnm_montecarlo_precision_lm(x, mu, sigma, mu_prior, sigma_prior, Binv, Z));
    return rcpp_result_gen;
END_RCPP
}
// c_obtain_moments_lrnm_montecarlo
Rcpp::List c_obtain_moments_lrnm_montecarlo(arma::mat Y, arma::vec mu, arma::mat sigma, arma::mat B, arma::mat& Z);
RcppExport SEXP _coda_count_c_obtain_moments_lrnm_montecarlo(SEXP YSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(c_obtain_moments_lrnm_montecarlo(Y, mu, sigma, B, Z));
    return rcpp_result_gen;
END_RCPP
}
// c_fit_lm_lrnm_montecarlo_centered
Rcpp::List c_fit_lm_lrnm_montecarlo_centered(arma::mat Y, arma::mat B, arma::mat X, arma::mat& Z, double eps, int max_iter);
RcppExport SEXP _coda_count_c_fit_lm_lrnm_montecarlo_centered(SEXP YSEXP, SEXP BSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_fit_lm_lrnm_montecarlo_centered(Y, B, X, Z, eps, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// ldnormal_vec
double ldnormal_vec(arma::vec h, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP _coda_count_ldnormal_vec(SEXP hSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ldnormal_vec(h, mu, inv_sigma));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_join_d1
arma::vec l_lrnm_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_join_d1(SEXP hSEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_join_d1(h, x, mu, inv_sigma, Binv));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_join_d2
arma::mat l_lrnm_join_d2(arma::vec h, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_join_d2(SEXP hSEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_join_d2(h, x, mu, inv_sigma, Binv));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_join_maximum
arma::vec l_lrnm_join_maximum(arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv, double eps, int max_iter);
RcppExport SEXP _coda_count_l_lrnm_join_maximum(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_join_maximum(x, mu, inv_sigma, Binv, eps, max_iter));
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
arma::mat c_rnormalSimplex(int n, arma::vec mu, arma::mat sigma, arma::mat Binv);
RcppExport SEXP _coda_count_c_rnormalSimplex(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rnormalSimplex(n, mu, sigma, Binv));
    return rcpp_result_gen;
END_RCPP
}
// c_rnormalmultinomial
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, arma::mat Binv);
RcppExport SEXP _coda_count_c_rnormalmultinomial(SEXP muSEXP, SEXP sigmaSEXP, SEXP sizeSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rnormalmultinomial(mu, sigma, size, Binv));
    return rcpp_result_gen;
END_RCPP
}
// pinv
arma::mat pinv(arma::mat X);
RcppExport SEXP _coda_count_pinv(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(pinv(X));
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
    {"_coda_count_H", (DL_FUNC) &_coda_count_H, 1},
    {"_coda_count_F", (DL_FUNC) &_coda_count_F, 1},
    {"_coda_count_c_ddm", (DL_FUNC) &_coda_count_c_ddm, 2},
    {"_coda_count_c_dm_fit", (DL_FUNC) &_coda_count_c_dm_fit, 3},
    {"_coda_count_c_dirichlet_alr_approx", (DL_FUNC) &_coda_count_c_dirichlet_alr_approx, 1},
    {"_coda_count_c_logistic_ilr_approximation", (DL_FUNC) &_coda_count_c_logistic_ilr_approximation, 1},
    {"_coda_count_c_logistic_alr_approximation", (DL_FUNC) &_coda_count_c_logistic_alr_approximation, 1},
    {"_coda_count_c_gaussian_product", (DL_FUNC) &_coda_count_c_gaussian_product, 2},
    {"_coda_count_c_gaussian_division", (DL_FUNC) &_coda_count_c_gaussian_division, 2},
    {"_coda_count_c_posterior_approximation", (DL_FUNC) &_coda_count_c_posterior_approximation, 4},
    {"_coda_count_c_fit_lrnm_gaussian_approx", (DL_FUNC) &_coda_count_c_fit_lrnm_gaussian_approx, 5},
    {"_coda_count_c_d_lrnm_hermite", (DL_FUNC) &_coda_count_c_d_lrnm_hermite, 5},
    {"_coda_count_c_moments_lrnm_hermite_precision_lm", (DL_FUNC) &_coda_count_c_moments_lrnm_hermite_precision_lm, 7},
    {"_coda_count_c_obtain_moments_lrnm_hermite", (DL_FUNC) &_coda_count_c_obtain_moments_lrnm_hermite, 5},
    {"_coda_count_c_fit_lm_lrnm_hermite_centered", (DL_FUNC) &_coda_count_c_fit_lm_lrnm_hermite_centered, 6},
    {"_coda_count_c_d_lrnm_montecarlo", (DL_FUNC) &_coda_count_c_d_lrnm_montecarlo, 5},
    {"_coda_count_c_moments_lrnm_montecarlo_precision_lm", (DL_FUNC) &_coda_count_c_moments_lrnm_montecarlo_precision_lm, 7},
    {"_coda_count_c_obtain_moments_lrnm_montecarlo", (DL_FUNC) &_coda_count_c_obtain_moments_lrnm_montecarlo, 5},
    {"_coda_count_c_fit_lm_lrnm_montecarlo_centered", (DL_FUNC) &_coda_count_c_fit_lm_lrnm_montecarlo_centered, 6},
    {"_coda_count_ldnormal_vec", (DL_FUNC) &_coda_count_ldnormal_vec, 3},
    {"_coda_count_l_lrnm_join_d1", (DL_FUNC) &_coda_count_l_lrnm_join_d1, 5},
    {"_coda_count_l_lrnm_join_d2", (DL_FUNC) &_coda_count_l_lrnm_join_d2, 5},
    {"_coda_count_l_lrnm_join_maximum", (DL_FUNC) &_coda_count_l_lrnm_join_maximum, 6},
    {"_coda_count_c_rmultinomial_Rcpp", (DL_FUNC) &_coda_count_c_rmultinomial_Rcpp, 2},
    {"_coda_count_c_rmultinomial", (DL_FUNC) &_coda_count_c_rmultinomial, 2},
    {"_coda_count_c_rdirichlet", (DL_FUNC) &_coda_count_c_rdirichlet, 2},
    {"_coda_count_c_rdirichletmultinomial", (DL_FUNC) &_coda_count_c_rdirichletmultinomial, 2},
    {"_coda_count_c_rnormal", (DL_FUNC) &_coda_count_c_rnormal, 3},
    {"_coda_count_c_rnormalSimplex", (DL_FUNC) &_coda_count_c_rnormalSimplex, 4},
    {"_coda_count_c_rnormalmultinomial", (DL_FUNC) &_coda_count_c_rnormalmultinomial, 4},
    {"_coda_count_pinv", (DL_FUNC) &_coda_count_pinv, 1},
    {"_coda_count_c_simplex_lattice", (DL_FUNC) &_coda_count_c_simplex_lattice, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_coda_count(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
