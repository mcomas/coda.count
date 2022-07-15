// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppEnsmallen.h>
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double c_d_lrnm_montecarlo(arma::vec x, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv,
                           arma::mat &Z){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_rows;

  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
  arma::mat N_posterior = c_posterior_approximation_vec(x, mu_prior, inv_sigma_prior, Binv);
  arma::vec mu = N_posterior.col(d);
  arma::mat sigma = N_posterior.head_cols(d);
  arma::mat inv_sigma = arma::inv_sympd(sigma);


  arma::mat Hz = Z * arma::chol(sigma);

  double l_cmult = l_multinomial_const(x);

  arma::vec h, p;
  double integral = 0;
  for(int i=0; i<n; i++){
    h = Hz.row(i).t() + mu;
    p = arma::exp(Binv * h);

    integral += exp(l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
      l_dnormal_vec(h, mu, inv_sigma) +
      l_multinomial(x, p/accu(p), l_cmult));
  }
  return integral/n;
}

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_montecarlo(arma::vec x,
                                    arma::vec mu, arma::mat sigma,
                                    arma::vec mu_prior, arma::mat inv_sigma_prior,
                                    arma::mat Binv, arma::mat &Z,
                                    arma::vec mu_centering){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_cols;

  arma::mat inv_sigma = arma::inv_sympd(sigma);

  arma::mat Hz = arma::chol(sigma).t() * Z;

  arma::vec h, p;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);

  arma::mat p_mu = arma::exp(Binv * mu);
  double constant = -arma::dot(log(p_mu/accu(p_mu)),x);
  for(int i=0; i<n; i++){
    h = Hz.col(i) + mu;
    p = arma::exp(Binv * h);
    double dens = exp(l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
                      l_dnormal_vec(h, mu, inv_sigma) +
                      arma::dot(log(p/accu(p)),x) + constant);
    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
  }

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_montecarlo_sigma_inverse(arma::vec x,
                                                  arma::vec mu, arma::mat inv_sigma,
                                                  arma::vec mu_prior, arma::mat inv_sigma_prior,
                                                  arma::mat Binv, arma::mat &Z,
                                                  arma::vec mu_centering){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_cols;

  arma::mat Hz = chol(arma::inv_sympd(inv_sigma)).t() * Z;

  arma::vec h, p;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);

  arma::mat p_mu = arma::exp(Binv * mu);
  double constant = -arma::dot(log(p_mu/accu(p_mu)),x);
  for(int i=0; i<n; i++){
    h = Hz.col(i) + mu;
    p = arma::exp(Binv * h);
    double dens = exp(l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
                      l_dnormal_vec(h, mu, inv_sigma) +
                      arma::dot(log(p/accu(p)),x) + constant);

    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
  }

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_cond_montecarlo(arma::vec x,
                                         arma::vec mu, arma::mat sigma,
                                         arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
                                         arma::mat Binv, arma::mat &Z,
                                         arma::vec mu_centering){

  unsigned dh1 =  Binv.n_cols-h2.size();

  arma::span I1 = arma::span(0, dh1-1);
  arma::span I2 = arma::span(dh1,  Binv.n_cols-1);

  arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));

  arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
  arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
  arma::mat inv_Sc = arma::inv_sympd(Sc);

  arma::mat inv_sigma = arma::inv_sympd(sigma);

  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
  double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);

  unsigned int index[dh1+1];
  for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
  int position = 0, k = 0;
  double M0 = 0;
  arma::vec M1 = arma::zeros(dh1);
  arma::mat M2 = arma::zeros(dh1,dh1);
  double l_cmult = l_multinomial_const(x);

  arma::mat Hz = arma::chol(sigma).t() * Z;
  unsigned n = Z.n_cols;
  // Rcpp::Rcout << "here";
  // unsigned d = x.n_elem - 1;
  //
  //
  // arma::mat inv_sigma = arma::inv_sympd(sigma);
  //
  // unsigned int index[d+1];
  // for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  // int position = 0, k = 0;
  //
  // arma::mat Hz = arma::chol(sigma).t() * Z;
  //
  arma::vec h1, h, p;
  // double M0 = 0;
  // arma::vec M1 = arma::zeros(d);
  // arma::mat M2 = arma::zeros(d,d);
  //
  arma::vec mu_c = join_cols(mu, h2);
  arma::mat p_mu = arma::exp(Binv * mu_c);
  double constant = -arma::dot(log(p_mu/accu(p_mu)),x);
  for(int i=0; i<n; i++){
    h1 = Hz.col(i) + mu;
    h = join_cols(h1, h2);
    p = arma::exp(Binv * h);
    double dens = exp(l_dnormal_vec(h1, Mc, inv_Sc) -
                      l_dnormal_vec(h1, mu, inv_sigma) +
                      arma::dot(log(p/accu(p)),x) + constant);
    M0 += dens;
    M1 += h1 * dens;
    M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
  }

  arma::mat moments(dh1, dh1+1);
  moments.col(dh1) = M1/M0;
  moments.head_cols(dh1) = M2/M0;
  return moments;
}

arma::mat c_moments_lrnm_montecarlo_old(arma::vec x,
                                        arma::vec mu, arma::mat sigma,
                                        arma::vec mu_prior, arma::mat inv_sigma_prior,
                                        arma::mat Binv, arma::mat &Z,
                                        arma::vec mu_centering){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_cols;

  arma::mat inv_sigma = arma::inv_sympd(sigma);

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;

  arma::mat Hz = arma::chol(sigma).t() * Z;

  // double l_cmult = l_multinomial_const(x);

  arma::vec h, p;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);

  arma::mat INV_SIGMA = inv_sigma_prior - inv_sigma;
  // arma::mat MU = inv(INV_SIGMA) * (inv_sigma_prior * mu_prior - inv_sigma * mu);
  // Rcpp::Rcout << INV_SIGMA;
  arma::mat MU = solve(INV_SIGMA, inv_sigma_prior * mu_prior - inv_sigma * mu);

  arma::mat p_mu = arma::exp(Binv * mu);
  // Rcpp::Rcout << p_mu;
  // Rcpp::Rcout << x;
  // if(FALSE & x[0] == 0 & x[1] == 12 & x[2] == 38){
  //   Rcpp::Rcout << "Start info 81" << std::endl;
  //   Rcpp::Rcout << "inv_sigma_prior:" << std::endl << inv_sigma_prior;
  //   Rcpp::Rcout << "inv_sigma:" << std::endl << inv_sigma;
  //   Rcpp::Rcout << p_mu << std::endl << MU << std::endl << INV_SIGMA << std::endl;
  //   Rcpp::Rcout << inv_sigma_prior * mu_prior - inv_sigma * mu;
  //   Rcpp::Rcout << "End info 81" << std::endl;
  // }
  double constant = -arma::dot(log(p_mu/accu(p_mu)),x);
  for(int i=0; i<n; i++){
    h = Hz.col(i) + mu;
    p = arma::exp(Binv * h);
    // Rcpp::Rcout << "(" << l_dnormal_prop_vec(h, MU, INV_SIGMA) << ", " << arma::dot(log(p/accu(p)),x) << ") = ";
    // Rcpp::Rcout << l_dnormal_prop_vec(h, MU, INV_SIGMA) + arma::dot(log(p/accu(p)),x) + constant << std::endl;
    double dens = exp(l_dnormal_vec(h, MU, INV_SIGMA) +
                      arma::dot(log(p/accu(p)),x) + constant);
    // double dens = exp(l_dnormal_prop_vec(h, mu_prior, inv_sigma_prior) -
    //                   l_dnormal_prop_vec(h, mu, inv_sigma) +
    //                   arma::dot(log(p/accu(p)),x));

    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
    // Rcpp::Rcout << M0 << M1 << M2;
  }
  // Rcpp::Rcout << M0 << M1 << M2 << std::endl;

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}


// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_obtain_moments_lrnm_montecarlo(arma::mat Y,
//                                             arma::vec mu, arma::mat sigma,
//                                             arma::mat B, arma::mat &Z){
//   int n = Y.n_rows;
//   int d = Y.n_cols - 1;
//
//   arma::mat Binv = pinv(B).t();
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//
//   arma::mat M1 = arma::zeros(d, n);
//   arma::cube M2 = arma::zeros(d, d, n);
//
//   for(int i = 0; i < Y.n_rows; i++){
//     arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu, inv_sigma, Binv);
//     arma::mat moments = c_moments_lrnm_montecarlo_precision_lm(Y.row(i).t(),
//                                                                N_posterior.col(d), N_posterior.head_cols(d),
//                                                                mu, inv_sigma,
//                                                                Binv, Z);
//     M1.col(i) = moments.col(d);
//     M2.slice(i) = moments.head_cols(d);
//   }
//   return Rcpp::List::create(M1, M2);
// }

//' @export
// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_lm_montecarlo(arma::mat Y, arma::mat B, arma::mat X, arma::mat &Z,
                                    double eps, int max_iter, arma::mat H0){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;

  arma::mat Zt = Z.t();


  arma::mat Binv = pinv(B).t();

  arma::mat H = H0;
  arma::mat invXtXX = arma::inv_sympd(X.t() * X) * X.t();
  arma::mat beta = invXtXX * H, beta_prev;
  arma::mat R = H - X * beta;
  arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;


  int current_iter = 0;
  arma::vec eigval;
  arma::mat eigvec;
  do{
    Rcpp::checkUserInterrupt();

    current_iter++;

    bool cor = eig_sym(eigval, eigvec, sigma_lm);

    if(eigval.min() < 1e-10 | eigval.max() < 1e-5){
      Rcpp::Rcout << "Covariance matrix is degenerate" << std::endl;
      return Rcpp::List::create(beta, sigma_lm, H, current_iter, eigval, sigma_lm);
    }

    inv_sigma_lm = arma::inv_sympd(sigma_lm);

    beta_prev = arma::mat(beta);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);

    for(int i = 0; i < Y.n_rows; i++){
      if(i % 50 == 0){
        Rcpp::checkUserInterrupt();
      }

      arma::mat mu_i =  X.row(i) * beta;

      arma::mat N_posterior_si = c_posterior_approximation_vec_sigma_inverse(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);

      arma::mat moments = c_moments_lrnm_montecarlo_sigma_inverse(Y.row(i).t(),
                                                                  N_posterior_si.col(d), N_posterior_si.head_cols(d),
                                                                  mu_i.t(), inv_sigma_lm,
                                                                  Binv, Zt, mu_i.t());


      H.row(i) = moments.col(d).t();
      M1 += moments.col(d);
      M2 += moments.head_cols(d);
    }
    beta = invXtXX * H;
    //R = H - X * beta;
    sigma_lm = M2 / n; //R.t() * R / (n-k);
    // Rcpp::Rcout << beta << norm(beta-beta_prev, 2) << " " << norm(beta-beta_prev, 1) << std::endl;
    //sigma = M2 / n - mu * mu.t();
  } while ( norm(beta-beta_prev, 1) > eps && current_iter < max_iter);
  // Last iteration
  inv_sigma_lm = arma::inv_sympd(sigma_lm);
  for(int i = 0; i < Y.n_rows; i++){
    arma::mat mu_i =  X.row(i) * beta;
    arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
    arma::mat moments = c_moments_lrnm_montecarlo(Y.row(i).t(),
                                                  N_posterior.col(d), N_posterior.head_cols(d),
                                                  mu_i.t(), inv_sigma_lm,
                                                  Binv, Zt, mu_i.t());
    H.row(i) = moments.col(d).t();
  }
  beta = invXtXX * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}
