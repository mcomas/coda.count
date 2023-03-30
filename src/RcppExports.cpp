// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEnsmallen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
Rcpp::List c_dm_fit(arma::mat& X, double eps, int maxiter);
RcppExport SEXP _coda_count_c_dm_fit(SEXP XSEXP, SEXP epsSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dm_fit(X, eps, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// hermite
arma::mat hermite(int order);
RcppExport SEXP _coda_count_hermite(SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(hermite(order));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_posterior_approximation_vec
arma::mat c_lrnm_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv, double eps, int niter);
RcppExport SEXP _coda_count_c_lrnm_posterior_approximation_vec(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP, SEXP epsSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_posterior_approximation_vec(x, mu, inv_sigma, Binv, eps, niter));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_posterior_approximation_vec_sigma_inverse
arma::mat c_lrnm_posterior_approximation_vec_sigma_inverse(arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv, double eps, int niter);
RcppExport SEXP _coda_count_c_lrnm_posterior_approximation_vec_sigma_inverse(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP, SEXP epsSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_posterior_approximation_vec_sigma_inverse(x, mu, inv_sigma, Binv, eps, niter));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_cond_posterior_approximation_vec
arma::mat c_lrnm_cond_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::vec h2, arma::mat& Binv, double eps, int niter);
RcppExport SEXP _coda_count_c_lrnm_cond_posterior_approximation_vec(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP h2SEXP, SEXP BinvSEXP, SEXP epsSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_cond_posterior_approximation_vec(x, mu, inv_sigma, h2, Binv, eps, niter));
    return rcpp_result_gen;
END_RCPP
}
// c_d_lrnm_hermite_mat
arma::vec c_d_lrnm_hermite_mat(arma::mat X, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_d_lrnm_hermite_mat(SEXP XSEXP, SEXP mu_priorSEXP, SEXP sigma_priorSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior(mu_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_hermite_mat(X, mu_prior, sigma_prior, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_posterior_moments_hermite
List c_lrnm_posterior_moments_hermite(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma, int order);
RcppExport SEXP _coda_count_c_lrnm_posterior_moments_hermite(SEXP XSEXP, SEXP clr_muSEXP, SEXP clr_sigmaSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clr_mu(clr_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type clr_sigma(clr_sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_posterior_moments_hermite(X, clr_mu, clr_sigma, order));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_hermite
List c_lrnm_fit_hermite(arma::mat& X, int order, double em_eps, int em_max_iter);
RcppExport SEXP _coda_count_c_lrnm_fit_hermite(SEXP XSEXP, SEXP orderSEXP, SEXP em_epsSEXP, SEXP em_max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< double >::type em_eps(em_epsSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_iter(em_max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_hermite(X, order, em_eps, em_max_iter));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_posterior_laplace_approximation
List c_lrnm_posterior_laplace_approximation(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma);
RcppExport SEXP _coda_count_c_lrnm_posterior_laplace_approximation(SEXP XSEXP, SEXP clr_muSEXP, SEXP clr_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clr_mu(clr_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type clr_sigma(clr_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_posterior_laplace_approximation(X, clr_mu, clr_sigma));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_laplace
List c_lrnm_fit_laplace(arma::mat& X, double em_eps, int em_max_iter);
RcppExport SEXP _coda_count_c_lrnm_fit_laplace(SEXP XSEXP, SEXP em_epsSEXP, SEXP em_max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type em_eps(em_epsSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_iter(em_max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_laplace(X, em_eps, em_max_iter));
    return rcpp_result_gen;
END_RCPP
}
// c_d_lrnm_montecarlo
arma::vec c_d_lrnm_montecarlo(arma::mat& X, arma::vec& mu, arma::mat& sigma, arma::mat& Binv, arma::mat& Z);
RcppExport SEXP _coda_count_c_d_lrnm_montecarlo(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BinvSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_montecarlo(X, mu, sigma, Binv, Z));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_posterior_moments_montecarlo
List c_lrnm_posterior_moments_montecarlo(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma, arma::mat& Z);
RcppExport SEXP _coda_count_c_lrnm_posterior_moments_montecarlo(SEXP XSEXP, SEXP clr_muSEXP, SEXP clr_sigmaSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clr_mu(clr_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type clr_sigma(clr_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_posterior_moments_montecarlo(X, clr_mu, clr_sigma, Z));
    return rcpp_result_gen;
END_RCPP
}
// c_lrnm_fit_montecarlo
List c_lrnm_fit_montecarlo(arma::mat& X, arma::mat& Z, double em_eps, int em_max_iter);
RcppExport SEXP _coda_count_c_lrnm_fit_montecarlo(SEXP XSEXP, SEXP ZSEXP, SEXP em_epsSEXP, SEXP em_max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type em_eps(em_epsSEXP);
    Rcpp::traits::input_parameter< int >::type em_max_iter(em_max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_lrnm_fit_montecarlo(X, Z, em_eps, em_max_iter));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_join_no_constant_vec
double l_lrnm_join_no_constant_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_join_no_constant_vec(SEXP hSEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_join_no_constant_vec(h, x, mu, inv_sigma, Binv));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_join_vec
double l_lrnm_join_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_join_vec(SEXP hSEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_join_vec(h, x, mu, inv_sigma, Binv));
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
// l_lrnm_cond_join_d1
arma::vec l_lrnm_cond_join_d1(arma::vec h1, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::vec h2, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_cond_join_d1(SEXP h1SEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP h2SEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_cond_join_d1(h1, x, mu, inv_sigma, h2, Binv));
    return rcpp_result_gen;
END_RCPP
}
// l_lrnm_cond_join_d2
arma::mat l_lrnm_cond_join_d2(arma::vec h1, arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::vec h2, arma::mat& Binv);
RcppExport SEXP _coda_count_l_lrnm_cond_join_d2(SEXP h1SEXP, SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP h2SEXP, SEXP BinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_cond_join_d2(h1, x, mu, inv_sigma, h2, Binv));
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
// l_lrnm_cond_join_maximum
arma::vec l_lrnm_cond_join_maximum(arma::vec x, arma::vec mu, arma::mat& inv_sigma, arma::vec h2, arma::mat& Binv, double eps, int max_iter);
RcppExport SEXP _coda_count_l_lrnm_cond_join_maximum(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP h2SEXP, SEXP BinvSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lrnm_cond_join_maximum(x, mu, inv_sigma, h2, Binv, eps, max_iter));
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
// pinv_sympd
arma::mat pinv_sympd(arma::mat X);
RcppExport SEXP _coda_count_pinv_sympd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(pinv_sympd(X));
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
    {"_coda_count_hermite", (DL_FUNC) &_coda_count_hermite, 1},
    {"_coda_count_c_lrnm_posterior_approximation_vec", (DL_FUNC) &_coda_count_c_lrnm_posterior_approximation_vec, 6},
    {"_coda_count_c_lrnm_posterior_approximation_vec_sigma_inverse", (DL_FUNC) &_coda_count_c_lrnm_posterior_approximation_vec_sigma_inverse, 6},
    {"_coda_count_c_lrnm_cond_posterior_approximation_vec", (DL_FUNC) &_coda_count_c_lrnm_cond_posterior_approximation_vec, 7},
    {"_coda_count_c_d_lrnm_hermite_mat", (DL_FUNC) &_coda_count_c_d_lrnm_hermite_mat, 5},
    {"_coda_count_c_lrnm_posterior_moments_hermite", (DL_FUNC) &_coda_count_c_lrnm_posterior_moments_hermite, 4},
    {"_coda_count_c_lrnm_fit_hermite", (DL_FUNC) &_coda_count_c_lrnm_fit_hermite, 4},
    {"_coda_count_c_lrnm_posterior_laplace_approximation", (DL_FUNC) &_coda_count_c_lrnm_posterior_laplace_approximation, 3},
    {"_coda_count_c_lrnm_fit_laplace", (DL_FUNC) &_coda_count_c_lrnm_fit_laplace, 3},
    {"_coda_count_c_d_lrnm_montecarlo", (DL_FUNC) &_coda_count_c_d_lrnm_montecarlo, 5},
    {"_coda_count_c_lrnm_posterior_moments_montecarlo", (DL_FUNC) &_coda_count_c_lrnm_posterior_moments_montecarlo, 4},
    {"_coda_count_c_lrnm_fit_montecarlo", (DL_FUNC) &_coda_count_c_lrnm_fit_montecarlo, 4},
    {"_coda_count_l_lrnm_join_no_constant_vec", (DL_FUNC) &_coda_count_l_lrnm_join_no_constant_vec, 5},
    {"_coda_count_l_lrnm_join_vec", (DL_FUNC) &_coda_count_l_lrnm_join_vec, 5},
    {"_coda_count_l_lrnm_join_d1", (DL_FUNC) &_coda_count_l_lrnm_join_d1, 5},
    {"_coda_count_l_lrnm_join_d2", (DL_FUNC) &_coda_count_l_lrnm_join_d2, 5},
    {"_coda_count_l_lrnm_cond_join_d1", (DL_FUNC) &_coda_count_l_lrnm_cond_join_d1, 6},
    {"_coda_count_l_lrnm_cond_join_d2", (DL_FUNC) &_coda_count_l_lrnm_cond_join_d2, 6},
    {"_coda_count_l_lrnm_join_maximum", (DL_FUNC) &_coda_count_l_lrnm_join_maximum, 6},
    {"_coda_count_l_lrnm_cond_join_maximum", (DL_FUNC) &_coda_count_l_lrnm_cond_join_maximum, 7},
    {"_coda_count_c_rmultinomial_Rcpp", (DL_FUNC) &_coda_count_c_rmultinomial_Rcpp, 2},
    {"_coda_count_c_rmultinomial", (DL_FUNC) &_coda_count_c_rmultinomial, 2},
    {"_coda_count_c_rdirichlet", (DL_FUNC) &_coda_count_c_rdirichlet, 2},
    {"_coda_count_c_rdirichletmultinomial", (DL_FUNC) &_coda_count_c_rdirichletmultinomial, 2},
    {"_coda_count_c_rnormal", (DL_FUNC) &_coda_count_c_rnormal, 3},
    {"_coda_count_c_rnormalSimplex", (DL_FUNC) &_coda_count_c_rnormalSimplex, 4},
    {"_coda_count_c_rnormalmultinomial", (DL_FUNC) &_coda_count_c_rnormalmultinomial, 4},
    {"_coda_count_pinv_sympd", (DL_FUNC) &_coda_count_pinv_sympd, 1},
    {"_coda_count_pinv", (DL_FUNC) &_coda_count_pinv, 1},
    {"_coda_count_c_simplex_lattice", (DL_FUNC) &_coda_count_c_simplex_lattice, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_coda_count(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
