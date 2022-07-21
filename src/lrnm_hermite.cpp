// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "hermite.h"
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// //' @export
// // [[Rcpp::export]]
// double c_d_lrnm_hermite(arma::vec x,
//                         arma::vec mu_prior, arma::mat sigma_prior,
//                         arma::mat Binv, int order){
//
//   unsigned d = Binv.n_cols;
//
//
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
//   arma::mat N_posterior = c_posterior_approximation_vec(x, mu_prior, inv_sigma_prior, Binv);
//   arma::vec mu = N_posterior.col(d);
//   arma::mat sigma = N_posterior.head_cols(d);
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   eig_sym(eigval, eigvec, sigma);
//   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
//
//   unsigned int index[d+1];
//   for(unsigned int i = 0; i <= d; i++) index[i] = 0;
//   int position = 0, k = 0;
//   double integral = 0;
//   double l_cmult = l_multinomial_const(x);
//
//   do{
//     double w = 0;
//     arma::vec h(d);
//     for(unsigned int i = 0; i < d; i++){
//       h(i) = uni_hermite(index[i],0);
//       w += uni_hermite(index[i],1);
//     }
//     h = mu + rotation * h;
//     arma::vec p = exp(Binv * h);
//     integral += exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
//       l_dnormal_vec(h, mu, inv_sigma) +
//       l_multinomial(x, p/accu(p), l_cmult));
//     // Calculate next coordinate
//     index[position]++;
//     while(index[position] == order){
//       index[position] = 0;
//       position++;
//       index[position]++;
//     }
//     position = 0;
//     k++;
//   } while (index[d] == 0);
//
//   return integral;
// }
//
// //' @export
// // [[Rcpp::export]]
// double c_d_lrnm_cond_hermite(arma::vec x,
//                              arma::vec mu_prior, arma::mat sigma_prior,
//                              arma::vec h2,
//                              arma::mat Binv, int order){
//
//  //unsigned d = Binv.n_cols;
//   unsigned dh1 =  Binv.n_cols-h2.size();
//
//   arma::span I1 = arma::span(0, dh1-1);
//   arma::span I2 = arma::span(dh1,  Binv.n_cols-1);
//
//   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
//
//   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
//   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   arma::mat inv_Sc = arma::inv_sympd(Sc);
//
//
//   arma::mat N_posterior = c_lrnm_cond_posterior_approximation_vec(x, Mc, inv_Sc, h2, Binv);
//
//   arma::vec mu = N_posterior.col(dh1);
//   arma::mat sigma = N_posterior.head_cols(dh1);
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   eig_sym(eigval, eigvec, sigma);
//   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
//
//   unsigned int index[dh1+1];
//   for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
//   int position = 0, k = 0;
//   double integral = 0;
//   double l_cmult = l_multinomial_const(x);
//
//   do{
//     double w = 0;
//     arma::vec h1(dh1);
//     for(unsigned int i = 0; i < dh1; i++){
//       h1(i) = uni_hermite(index[i],0);
//       w += uni_hermite(index[i],1);
//     }
//     h1 = mu + rotation * h1;
//     arma::vec h = join_cols(h1, h2);
//
//     arma::vec p = exp(Binv * h);
//     integral += exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
//       l_dnormal_vec(h1, mu, inv_sigma) +
//       l_multinomial(x, p/accu(p), l_cmult));
//
//     // Calculate next coordinate
//     index[position]++;
//     while(index[position] == order){
//       index[position] = 0;
//       position++;
//       index[position]++;
//     }
//     position = 0;
//     k++;
//   } while (index[dh1] == 0);
//
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
//   return integral * exp(l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior));
// }
//
//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_hermite(arma::vec x,
                                 arma::vec mu, arma::mat sigma,
                                 arma::vec mu_prior, arma::mat sigma_prior,
                                 arma::mat Binv, int order,
                                 arma::vec mu_centering){
  unsigned d = Binv.n_cols;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  arma::mat inv_sigma = arma::inv_sympd(sigma);
  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);
  double l_cmult = l_multinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    arma::vec p = exp(Binv * h);
    double dens = exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
                      l_dnormal_vec(h, mu, inv_sigma) +
                      l_multinomial(x, p/accu(p), l_cmult));
    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_centering) * (h-mu_centering).t() * dens;

    // Calculate next coordinate
    index[position]++;
    while(index[position] == order){
      index[position] = 0;
      position++;
      index[position]++;
    }
    position = 0;
    k++;
  } while (index[d] == 0);

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}

// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_hermite_sigma_inverse(arma::vec x,
//                                                arma::vec mu, arma::mat inv_sigma,
//                                                arma::vec mu_prior, arma::mat inv_sigma_prior,
//                                                arma::mat Binv, int order,
//                                                arma::vec mu_centering){
//   unsigned d = Binv.n_cols;
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   eig_sym(eigval, eigvec, inv_sigma);
//   arma::mat rotation = eigvec * arma::diagmat(sqrt(1/eigval));
//
//   unsigned int index[d+1];
//   for(unsigned int i = 0; i <= d; i++) index[i] = 0;
//   int position = 0, k = 0;
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(d);
//   arma::mat M2 = arma::zeros(d,d);
//   double l_cmult = l_multinomial_const(x);
//
//   do{
//     double w = 0;
//     arma::vec h(d);
//     for(unsigned int i = 0; i < d; i++){
//       h(i) = uni_hermite(index[i],0);
//       w += uni_hermite(index[i],1);
//     }
//     h = mu + rotation * h;
//     arma::vec p = exp(Binv * h);
//     double dens = exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
//                       l_dnormal_vec(h, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult));
//     M0 += dens;
//     M1 += h * dens;
//     M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
//
//     // Calculate next coordinate
//     index[position]++;
//     while(index[position] == order){
//       index[position] = 0;
//       position++;
//       index[position]++;
//     }
//     position = 0;
//     k++;
//   } while (index[d] == 0);
//
//   arma::mat moments(d, d+1);
//   moments.col(d) = M1/M0;
//   moments.head_cols(d) = M2/M0;
//   return moments;
// }

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_cond_hermite_1d(arma::vec x,
                                         arma::vec mu, double sigma,
                                         arma::vec Mc, arma::mat inv_Sc, arma::vec h2,
                                         arma::mat Binv, int order,
                                         double mu_centering){


  arma::span I1 = arma::span(0, 0);
  arma::span I2 = arma::span(1,  Binv.n_cols-1);

  // arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);
  //
  // arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
  // arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
  // arma::mat inv_Sc = arma::inv_sympd(Sc, arma::inv_opts::allow_approx);

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  arma::mat inv_sigma(1,1);
  inv_sigma(0,0) = 1/sigma;

  int position = 0, k = 0;
  double M0 = 0;
  arma::vec M1 = arma::zeros(1);
  arma::mat M2 = arma::zeros(1, 1);
  double l_cmult = l_multinomial_const(x);
  arma::vec h_mu = arma::vec(1+h2.size());
  h_mu(I1) = mu;
  h_mu(I2) = h2;
  arma::vec p_mu = exp(Binv * h_mu);

  double w;
  arma::vec h1(1);
  arma::vec h = arma::vec(1+h2.size());
  h(I2) = h2;
  for(int i = 0; i < order; i++){

    h1 = mu + sqrt(sigma) * uni_hermite(i,0);
    w  = uni_hermite(i,1);

    h(I1) = h1;
    arma::vec p = exp(Binv * h);
    // if(false){
    //   Rcpp::Rcout << "l_dnormal_vec(h1, Mc, inv_Sc):" << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
    //   Rcpp::Rcout << "l_dnormal_vec(h1, mu, inv_sigma):" << l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
    //   Rcpp::Rcout << "l_multinomial(x, p/accu(p), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
    //   Rcpp::Rcout << "l_multinomial(x, p_mu/accu(p_mu), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
    //   Rcpp::Rcout << "ldens: " << w + l_dnormal_vec(h1, Mc, inv_Sc) -
    //     l_dnormal_vec(h1, mu, inv_sigma) +
    //     l_multinomial(x, p/accu(p), l_cmult) -
    //     l_multinomial(x, p_mu/accu(p_mu), l_cmult) << std::endl << std::endl;
    // }
    double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
                      l_dnormal_vec(h1, mu, inv_sigma) +
                      l_multinomial(x, p/accu(p), l_cmult) -
                      l_multinomial(x, p_mu/accu(p_mu), l_cmult));  // For big x, low variability.


    M0 += dens;
    M1 += h1 * dens;
    M2 += (h1-mu_centering) * (h1-mu_centering) * dens;
    // Calculate next coordinate
  }


  arma::mat moments(1, 2);
  moments.col(1) = M1/M0;
  moments.head_cols(1) = M2/M0;
  return moments;
}

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_cond_hermite(arma::vec x,
                                      arma::vec mu, arma::mat sigma,
                                      arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
                                      arma::mat Binv, int order,
                                      arma::vec mu_centering){
  unsigned dh1 =  Binv.n_cols-h2.size();

  arma::span I1 = arma::span(0, dh1-1);
  arma::span I2 = arma::span(dh1,  Binv.n_cols-1);

  arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);

  arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
  arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
  arma::mat inv_Sc = arma::inv_sympd(Sc, arma::inv_opts::allow_approx);

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  arma::mat inv_sigma = arma::inv_sympd(sigma, arma::inv_opts::allow_approx);

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  Rcpp::Rcout << rotation << std::endl;
  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);
  double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);

  unsigned int index[dh1+1];
  for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
  int position = 0, k = 0;
  double M0 = 0;
  arma::vec M1 = arma::zeros(dh1);
  arma::mat M2 = arma::zeros(dh1,dh1);
  double l_cmult = l_multinomial_const(x);
  arma::vec h_mu = join_cols(mu, h2);
  arma::vec p_mu = exp(Binv * h_mu);
  do{
    double w = 0;
    arma::vec h1(dh1);
    for(unsigned int i = 0; i < dh1; i++){
      h1(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }

    h1 = mu + rotation * h1; // rotation 0x0

    arma::vec h = join_cols(h1, h2);
    arma::vec p = exp(Binv * h);
    // Rcpp::Rcout << "l_ph2:" << l_ph2 << std::endl;
    // Rcpp::Rcout << "l_dnormal_vec(h1, Mc, inv_Sc):" << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
    // Rcpp::Rcout << "l_dnormal_vec(h1, mu, inv_sigma):" << l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
    // Rcpp::Rcout << "l_multinomial(x, p/accu(p), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
    // Rcpp::Rcout << "l_multinomial(x, p_mu/accu(p_mu), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
    double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
                      l_dnormal_vec(h1, mu, inv_sigma) +
                      l_multinomial(x, p/accu(p), l_cmult) -
                      l_multinomial(x, p_mu/accu(p_mu), l_cmult));  // For big x, low variability.

    // Rcpp::Rcout << "ldens: " << w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
    //   l_dnormal_vec(h1, mu, inv_sigma) +
    //   l_multinomial(x, p/accu(p), l_cmult) << std::endl << std::endl;
    M0 += dens;
    M1 += h1 * dens;
    M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
    // Calculate next coordinate
    index[position]++;
    while(index[position] == order){
      index[position] = 0;
      position++;
      index[position]++;
    }
    position = 0;
    k++;
  } while (index[dh1] == 0);


  arma::mat moments(dh1, dh1+1);
  moments.col(dh1) = M1/M0;
  moments.head_cols(dh1) = M2/M0;
  return moments;
}
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_cond_hermite2(arma::vec x,
//                                        arma::vec mu, arma::mat sigma,
//                                        arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
//                                        arma::mat Binv, int order,
//                                        arma::vec mu_centering){
//   unsigned dh1 =  Binv.n_cols-h2.size();
//
//   arma::span I1 = arma::span(0, dh1-1);
//   arma::span I2 = arma::span(dh1,  Binv.n_cols-1);
//
//   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
//
//   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
//   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   arma::mat inv_Sc = arma::inv_sympd(Sc);
//
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//
//   eig_sym(eigval, eigvec, sigma);
//   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
//
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
//   // double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);
//
//   unsigned int index[dh1+1];
//   for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
//   int position = 0, k = 0;
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(dh1);
//   arma::mat M2 = arma::zeros(dh1,dh1);
//   double l_cmult = l_multinomial_const(x);
//
//   do{
//     double w = 0;
//     arma::vec h1(dh1);
//     for(unsigned int i = 0; i < dh1; i++){
//       h1(i) = uni_hermite(index[i],0);
//       w += uni_hermite(index[i],1);
//     }
//     h1 = mu + rotation * h1;
//     arma::vec h = join_cols(h1, h2);
//     arma::vec p = exp(Binv * h);
//     double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
//                       l_dnormal_vec(h1, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult));
//     // Rcpp::Rcout << dens;
//     M0 += dens;
//     M1 += h1 * dens;
//     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
//     // Calculate next coordinate
//     index[position]++;
//     while(index[position] == order){
//       index[position] = 0;
//       position++;
//       index[position]++;
//     }
//     position = 0;
//     k++;
//   } while (index[dh1] == 0);
//
//   arma::mat moments(dh1, dh1+1);
//   moments.col(dh1) = M1/M0;
//   moments.head_cols(dh1) = M2/M0;
//   return moments;
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_cond1_hermite(arma::vec h, int i1,
//                                        arma::vec x, arma::vec mu, double sigma,
//                                        arma::vec mu_prior, arma::mat sigma_prior,
//                                        arma::mat Binv, int order,
//                                        arma::vec mu_centering){
//
//
//   arma::uvec I1(1);
//   I1(0) = i1 - 1;
//   arma::uvec I2(Binv.n_cols-1);
//   for(int i = 0, k = 0; i < Binv.n_cols; i++){
//     if(i + 1 != i1){
//       I2(k) = i;
//       k++;
//     }
//   }
//   // Rcpp::Rcout << I1 << std::endl;
//   // Rcpp::Rcout << I2 << std::endl;
//
//   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
//
//   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h(I2)-mu_prior(I2));
//   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   arma::mat inv_Sc = arma::inv_sympd(Sc);
//
//
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::mat inv_sigma(1,1);
//   inv_sigma(0,0) = 1/sigma;
//
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
//   double l_ph2 = l_dnormal_vec(h(I2), mu_prior(I2), inv_sigma_prior_2);
//
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(1);
//   arma::mat M2 = arma::zeros(1,1);
//   double l_cmult = l_multinomial_const(x);
//
//   // Rcpp::Rcout << uni_hermite << std::endl;
//   double w = 0;
//   arma::vec h1;
// //   arma::vec h = arma::zeros(Binv.n_cols);
// //   h(I2) = h2;
//   for(int i = 0; i < order; i++){
//     h1 = sqrt(sigma) * uni_hermite(i,0) + mu;
//     w = uni_hermite(i,1);
//     h(I1) = h1;
//     arma::vec p = exp(Binv * h);
//
//     double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
//                       l_dnormal_vec(h1, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult));
//
//     M0 += dens;
//     M1 += h1 * dens;
//     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
//   }
//
//
//
//   arma::mat moments(1, 2);
//   moments.col(1) = M1/M0;
//   moments.col(0) = M2/M0;
//   return moments;
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_indep_lrnm_cond_hermite(arma::vec x, int i1,
//                                             arma::vec mu, double sigma,
//                                             arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
//                                             arma::mat Binv, int order,
//                                             arma::vec mu_centering){
//   unsigned dh1 =  Binv.n_cols-h2.size();
//
//   arma::uvec I1(1);
//   I1(0) = i1 - 1;
//   arma::uvec I2(Binv.n_cols-1);
//
//   for(int i = 0, k = 0; i < Binv.n_cols; i++){
//     if(i + 1 != i1){
//       I2(k) = i;
//       k++;
//     }
//   }
//
//
//   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
//
//   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
//   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   arma::mat inv_Sc = arma::inv_sympd(Sc);
//
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   arma::mat inv_sigma(1,1);
//   inv_sigma(0,0) = 1/sigma;
//
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
//   double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);
//
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(dh1);
//   arma::mat M2 = arma::zeros(dh1,dh1);
//   double l_cmult = l_multinomial_const(x);
//
//   //Rcpp::Rcout << uni_hermite << std::endl;
//   double w = 0;
//   arma::vec h1;
//   arma::vec h = arma::zeros(Binv.n_cols);
//   h(I2) = h2;
//   for(int i = 0; i < order; i++){
//     h1 = sqrt(sigma) * uni_hermite(i,0) + mu;
//     w = uni_hermite(i,1);
//     h(I1) = h1;
//     arma::vec p = exp(Binv * h);
//
//     double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
//                       l_dnormal_vec(h1, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult));
//
//     M0 += dens;
//     M1 += h1 * dens;
//     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
//   }
//
//   /*
//    do{
//    double w = 0;
//    arma::vec h1(dh1);
//    for(unsigned int i = 0; i < dh1; i++){
//    h1(i) = uni_hermite(index[i],0);
//    w += uni_hermite(index[i],1);
//    }
//    h1 = mu + h1;
//    arma::vec h = join_cols(h1, h2);
//    arma::vec p = exp(Binv * h);
//    // Rcpp::Rcout << "Nh2: " << l_ph2 << std::endl;
//    // Rcpp::Rcout << "Nc: " << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
//    // Rcpp::Rcout << "N: " << -l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
//    // Rcpp::Rcout << "Mult: " << l_multinomial(x, p/accu(p), l_cmult) << std::endl << std::endl;
//    double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
//    l_dnormal_vec(h1, mu, inv_sigma) +
//    l_multinomial(x, p/accu(p), l_cmult));
//    // Rcpp::Rcout << dens;
//    M0 += dens;
//    M1 += h1 * dens;
//    M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
//    // Calculate next coordinate
//    index[position]++;
//    while(index[position] == order){
//    index[position] = 0;
//    position++;
//    index[position]++;
//    }
//    position = 0;
//    k++;
//    } while (index[dh1] == 0);
//    */
//
//   arma::mat moments(dh1, dh1+1);
//   moments.col(dh1) = M1/M0;
//   moments.head_cols(dh1) = M2/M0;
//   return moments;
// }
//
// // //' @export
// // // [[Rcpp::export]]
// // Rcpp::List c_obtain_moments_lrnm_hermite(arma::mat Y,
// //                                          arma::vec mu, arma::mat sigma,
// //                                          arma::mat B, int order){
// //   int n = Y.n_rows;
// //   int d = Y.n_cols - 1;
// //
// //   arma::mat Binv = pinv(B).t();
// //   arma::mat inv_sigma = arma::inv_sympd(sigma);
// //
// //   arma::mat M1 = arma::zeros(d, n);
// //   arma::cube M2 = arma::zeros(d, d, n);
// //
// //   for(int i = 0; i < Y.n_rows; i++){
// //     arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu, inv_sigma, Binv);
// //     arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
// //                                                             N_posterior.col(d), N_posterior.head_cols(d),
// //                                                             mu, sigma,
// //                                                             Binv, order);
// //     M1.col(i) = moments.col(d);
// //     M2.slice(i) = moments.head_cols(d);
// //   }
// //   return Rcpp::List::create(M1, M2);
// // }
//
// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_fit_lrnm_lm_hermite(arma::mat Y, arma::mat B, arma::mat X, int order,
//                                  double eps, int max_iter, arma::mat H0){
//
//   int n = Y.n_rows;
//   int k = X.n_cols;
//   int d = Y.n_cols - 1;
//
//   arma::mat Binv = pinv(B).t();
//   arma::mat H = H0;
//
//   arma::mat beta = arma::inv(X.t() * X) * X.t() * H, beta_prev;
//   arma::mat R = H - X * beta;
//   arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;
//
//
//   int current_iter = 0;
//   arma::vec eigval;
//   arma::mat eigvec;
//   do{
//     Rcpp::checkUserInterrupt();
//
//     current_iter++;
//
//     bool cor = eig_sym(eigval, eigvec, sigma_lm);
//     if(eigval.min() < 1e-10){
//       Rcpp::Rcout << "Covariance matrix is degenerate" << std::endl;
//       return Rcpp::List::create(beta, sigma_lm, H, current_iter, eigval, sigma_lm);
//     }
//     inv_sigma_lm = arma::inv_sympd(sigma_lm);
//     beta_prev = arma::mat(beta);
//     arma::vec M1 = arma::zeros(d);
//     arma::mat M2 = arma::zeros(d, d);
//     for(int i = 0; i < Y.n_rows; i++){
//       arma::mat mu_i =  X.row(i) * beta;
//       arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
//       arma::mat moments = c_moments_lrnm_hermite(Y.row(i).t(),
//                                                  N_posterior.col(d), N_posterior.head_cols(d),
//                                                  mu_i.t(), sigma_lm,
//                                                  Binv, order, mu_i.t());
//       H.row(i) = moments.col(d).t();
//       M1 += moments.col(d);
//       M2 += moments.head_cols(d);
//     }
//     beta = arma::inv(X.t() * X) * X.t() * H;
//     //R = H - X * beta;
//     sigma_lm = M2 / n; //R.t() * R / (n-k);
//     //sigma = M2 / n - mu * mu.t();
//   } while ( norm(beta-beta_prev, 1) > eps && current_iter < max_iter);
//
//   // Last iteration
//   inv_sigma_lm = arma::inv_sympd(sigma_lm);
//   for(int i = 0; i < Y.n_rows; i++){
//     arma::mat mu_i =  X.row(i) * beta;
//     arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
//     arma::mat moments = c_moments_lrnm_hermite(Y.row(i).t(),
//                                                N_posterior.col(d), N_posterior.head_cols(d),
//                                                mu_i.t(), sigma_lm,
//                                                Binv, order, mu_i.t());
//
//     H.row(i) = moments.col(d).t();
//   }
//   beta = arma::inv(X.t() * X) * X.t() * H;
//
//   return Rcpp::List::create(beta, sigma_lm, H, current_iter);
// }
