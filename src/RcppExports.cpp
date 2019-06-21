// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// c_d_lrnm_hermite
double c_d_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_d_lrnm_hermite(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_hermite(x, mu, sigma, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_m1_lrnm_hermite
arma::vec c_m1_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_m1_lrnm_hermite(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_m1_lrnm_hermite(x, mu, sigma, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_m2_lrnm_hermite
arma::mat c_m2_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv, int order);
RcppExport SEXP _coda_count_c_m2_lrnm_hermite(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BinvSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Binv(BinvSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(c_m2_lrnm_hermite(x, mu, sigma, Binv, order));
    return rcpp_result_gen;
END_RCPP
}
// c_d_lrnm_hermite_
double c_d_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma, int order, int step_by, double eps, int max_steps);
RcppExport SEXP _coda_count_c_d_lrnm_hermite_(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(c_d_lrnm_hermite_(x, mu, sigma, order, step_by, eps, max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_m1_lrnm_hermite_
arma::vec c_m1_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma, int order, int step_by, double eps, int max_steps);
RcppExport SEXP _coda_count_c_m1_lrnm_hermite_(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(c_m1_lrnm_hermite_(x, mu, sigma, order, step_by, eps, max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_m2_lrnm_hermite_
arma::mat c_m2_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma, int order, int step_by, double eps, int max_steps);
RcppExport SEXP _coda_count_c_m2_lrnm_hermite_(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(c_m2_lrnm_hermite_(x, mu, sigma, order, step_by, eps, max_steps));
    return rcpp_result_gen;
END_RCPP
}
// c_fit_lrnm_hermite_
Rcpp::List c_fit_lrnm_hermite_(arma::mat X, arma::vec mu0, arma::mat sigma0, int order, int step_by, double eps, int max_steps, int em_max_steps);
RcppExport SEXP _coda_count_c_fit_lrnm_hermite_(SEXP XSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP orderSEXP, SEXP step_bySEXP, SEXP epsSEXP, SEXP max_stepsSEXP, SEXP em_max_stepsSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(c_fit_lrnm_hermite_(X, mu0, sigma0, order, step_by, eps, max_steps, em_max_steps));
    return rcpp_result_gen;
END_RCPP
}
// ldnormal2
arma::vec ldnormal2(arma::mat H, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP _coda_count_ldnormal2(SEXP HSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ldnormal2(H, mu, inv_sigma));
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
// lpnm_join
arma::mat lpnm_join(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::mat P, arma::mat H);
RcppExport SEXP _coda_count_lpnm_join(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP PSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(lpnm_join(x, mu, inv_sigma, P, H));
    return rcpp_result_gen;
END_RCPP
}
// lpnm_join_no_constant
double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::vec p, arma::vec h);
RcppExport SEXP _coda_count_lpnm_join_no_constant(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP pSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(lpnm_join_no_constant(x, mu, inv_sigma, p, h));
    return rcpp_result_gen;
END_RCPP
}
// lpnm_join_maximum_alr
arma::vec lpnm_join_maximum_alr(arma::vec x, arma::vec mu_alr, arma::mat& inv_sigma_alr, arma::vec a, double eps, int max_iter);
RcppExport SEXP _coda_count_lpnm_join_maximum_alr(SEXP xSEXP, SEXP mu_alrSEXP, SEXP inv_sigma_alrSEXP, SEXP aSEXP, SEXP epsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_alr(mu_alrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inv_sigma_alr(inv_sigma_alrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(lpnm_join_maximum_alr(x, mu_alr, inv_sigma_alr, a, eps, max_iter));
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
    {"_coda_count_c_ddm", (DL_FUNC) &_coda_count_c_ddm, 2},
    {"_coda_count_c_dm_fit", (DL_FUNC) &_coda_count_c_dm_fit, 3},
    {"_coda_count_c_d_lrnm_hermite", (DL_FUNC) &_coda_count_c_d_lrnm_hermite, 5},
    {"_coda_count_c_m1_lrnm_hermite", (DL_FUNC) &_coda_count_c_m1_lrnm_hermite, 5},
    {"_coda_count_c_m2_lrnm_hermite", (DL_FUNC) &_coda_count_c_m2_lrnm_hermite, 5},
    {"_coda_count_c_d_lrnm_hermite_", (DL_FUNC) &_coda_count_c_d_lrnm_hermite_, 7},
    {"_coda_count_c_m1_lrnm_hermite_", (DL_FUNC) &_coda_count_c_m1_lrnm_hermite_, 7},
    {"_coda_count_c_m2_lrnm_hermite_", (DL_FUNC) &_coda_count_c_m2_lrnm_hermite_, 7},
    {"_coda_count_c_fit_lrnm_hermite_", (DL_FUNC) &_coda_count_c_fit_lrnm_hermite_, 8},
    {"_coda_count_ldnormal2", (DL_FUNC) &_coda_count_ldnormal2, 3},
    {"_coda_count_ldnormal", (DL_FUNC) &_coda_count_ldnormal, 3},
    {"_coda_count_lpmultinomial_const", (DL_FUNC) &_coda_count_lpmultinomial_const, 1},
    {"_coda_count_lpnm_join", (DL_FUNC) &_coda_count_lpnm_join, 5},
    {"_coda_count_lpnm_join_no_constant", (DL_FUNC) &_coda_count_lpnm_join_no_constant, 5},
    {"_coda_count_lpnm_join_maximum_alr", (DL_FUNC) &_coda_count_lpnm_join_maximum_alr, 6},
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
