// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "hermite.h"
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// [[Rcpp::export]]
double c_d_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                        int order){

  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  double integral = 0;
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
    integral += exp(w + l_multinomial(x, p/accu(p), l_cmult));
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

  return integral;
}


// [[Rcpp::export]]
arma::vec c_m1_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                            arma::mat Binv, int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  arma::vec integral = arma::zeros(d);
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
    integral += h * exp(w + l_multinomial(x, p/accu(p), l_cmult));

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

  return integral;
}


// [[Rcpp::export]]
arma::mat c_m2_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                            int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  arma::mat integral = arma::zeros(d,d);
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
    integral += h * h.t()* exp(w + l_multinomial(x, p/accu(p), l_cmult));

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

  return integral;
}

// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
//                                  int order){
//   unsigned d = x.n_elem - 1;
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
//     double dens = exp(w + l_multinomial(x, p/accu(p), l_cmult));
//     M0 += dens;
//     M1 += h * dens;
//     M2 += h * h.t() * dens;
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

// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_hermite_precision(arma::vec x,
//                                            arma::vec mu, arma::mat sigma,
//                                            arma::vec mu_prior, arma::mat sigma_prior,
//                                            arma::mat Binv, int order){
//   unsigned d = x.n_elem - 1;
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
//
//   eig_sym(eigval, eigvec, sigma);
//   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
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
//     double dens = exp(w + ldnormal_vec(h, mu_prior, inv_sigma_prior) -
//                       ldnormal_vec(h, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult));
//     M0 += dens;
//     M1 += h * dens;
//     M2 += h * h.t() * dens;
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


// [[Rcpp::export]]
arma::mat c_moments_lrnm_hermite_precision_lm(arma::vec x,
                                              arma::vec mu, arma::mat sigma,
                                              arma::vec mu_prior, arma::mat sigma_prior,
                                              arma::mat Binv, int order){
  unsigned d = x.n_elem - 1;
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
    double dens = exp(w + ldnormal_vec(h, mu_prior, inv_sigma_prior) -
                      ldnormal_vec(h, mu, inv_sigma) +
                      l_multinomial(x, p/accu(p), l_cmult));
    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_prior) * (h-mu_prior).t() * dens;

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
// Rcpp::List c_fit_lrnm_hermite_precision(arma::mat X, arma::mat B, int order,
//                                         double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){
//
//   int n = X.n_rows;
//   int d = X.n_cols - 1;
//   arma::vec alpha = c_dm_fit_alpha(X);
//   arma::mat Binv = pinv(B).t();
//   arma::mat P = arma::mat(X);
//   P.each_row() += alpha.t();
//
//   arma::mat H = arma::log(P) * B;
//   arma::vec mu = mean(H,0).t();
//   arma::vec mu_prev;
//   arma::mat sigma = cov(H);
//   int current_iter = 0;
//   do{
//     current_iter++;
//     arma::mat inv_sigma = arma::inv_sympd(sigma);
//     //Rcpp::Rcout << "Current: " << current_iter << std::endl;
//     mu_prev = arma::vec(mu);
//     arma::vec M1 = arma::zeros(d);
//     arma::mat M2 = arma::zeros(d, d);
//     for(int i = 0; i < H.n_rows; i++){
//       // Rcpp::Rcout << mu << std::endl;
//       // Rcpp::Rcout << sigma << std::endl;
//       arma::mat N_posterior = c_posterior_approximation(X.row(i).t(), mu, inv_sigma, Binv);
//       //Rcpp::Rcout << N_posterior;
//       arma::mat moments = c_moments_lrnm_hermite_precision(X.row(i).t(),
//                                                            N_posterior.col(d), N_posterior.head_cols(d),
//                                                            mu, sigma,
//                                                            Binv, order);
//       // Rcpp::Rcout << i << std::endl;
//       // Rcpp::Rcout << moments << std::endl;
//       M1 += moments.col(d);
//       M2 += moments.head_cols(d);
//       //Rcpp::Rcout << M1 << std::endl;
//       // Rcpp::Rcout << M2 << std::endl;
//     }
//     mu = M1 / n;
//     sigma = M2 / n - mu * mu.t();
//   } while ( norm(mu-mu_prev, 2) > eps && current_iter < max_iter);
//
//   return Rcpp::List::create(mu, sigma, current_iter);
// }

// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_fit_lrnm_hermite(arma::mat X, arma::mat B, int order,
//                               double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){
//
//   int n = X.n_rows;
//   int d = X.n_cols - 1;
//   arma::vec alpha = c_dm_fit_alpha(X);
//   arma::mat Binv = pinv(B).t();
//   arma::mat P = arma::mat(X);
//   P.each_row() += alpha.t();
//
//   arma::mat H = arma::log(P) * B;
//   arma::vec mu = mean(H,0).t();
//   arma::vec mu_prev;
//   arma::mat sigma = cov(H);
//   int current_iter = 0;
//   do{
//     current_iter++;
//     //Rcpp::Rcout << "Current: " << current_iter << std::endl;
//     mu_prev = arma::vec(mu);
//     arma::vec M1 = arma::zeros(d);
//     arma::mat M2 = arma::zeros(d, d);
//     for(int i = 0; i < H.n_rows; i++){
//
//       arma::mat moments = c_moments_lrnm_hermite(X.row(i).t(),
//                                                  mu, sigma,
//                                                  Binv, order);
//       // Rcpp::Rcout << i << std::endl;
//       // Rcpp::Rcout << moments << std::endl;
//       M1 += moments.col(d);
//       M2 += moments.head_cols(d);
//       //Rcpp::Rcout << M1 << std::endl;
//       // Rcpp::Rcout << M2 << std::endl;
//     }
//     mu = M1 / n;
//     sigma = M2 / n - mu * mu.t();
//   } while ( norm(mu-mu_prev, 2) > eps && current_iter < max_iter);
//
//   return Rcpp::List::create(mu, sigma, current_iter);
// }

// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_fit_lm_lrnm_hermite(arma::mat Y, arma::mat B, arma::mat X, int order,
//                                  double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){
//
//   int n = Y.n_rows;
//   int k = X.n_cols;
//   int d = Y.n_cols - 1;
//   arma::vec alpha = c_dm_fit_alpha(Y);
//   arma::mat Binv = pinv(B).t();
//   arma::mat P = arma::mat(Y);
//   P.each_row() += alpha.t();
//
//   arma::mat H = arma::log(P) * B;
//
//   arma::mat beta = arma::inv(X.t() * X) * X.t() * H, beta_prev;
//   arma::mat R = H - X * beta;
//   arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;
//
//
//   int current_iter = 0;
//   do{
//     current_iter++;
//     //Rcpp::Rcout << "sigma:" << std::endl << sigma_lm;
//     inv_sigma_lm = arma::inv_sympd(sigma_lm);
//     //Rcpp::Rcout << "inv sigma:" << std::endl << inv_sigma_lm;
//     beta_prev = arma::mat(beta);
//     arma::vec M1 = arma::zeros(d);
//     arma::mat M2 = arma::zeros(d, d);
//     for(int i = 0; i < Y.n_rows; i++){
//       arma::mat mu_i =  X.row(i) * beta;
//       //Rcpp::Rcout << "start optimisation";
//       //Rcpp::Rcout << mu_i << inv_sigma_lm;
//       arma::mat N_posterior = c_posterior_approximation(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
//       arma::mat moments = c_moments_lrnm_hermite_precision(Y.row(i).t(),
//                                                            N_posterior.col(d), N_posterior.head_cols(d),
//                                                            mu_i.t(), sigma_lm,
//                                                            Binv, order);
//       H.row(i) = moments.col(d).t();
//       M1 += moments.col(d);
//       M2 += moments.head_cols(d);
//     }
//     //Rcpp::Rcout << "ending" << std::endl;
//     beta = arma::inv(X.t() * X) * X.t() * H;
//     //R = H - X * beta;
//     sigma_lm = M2 / n; //R.t() * R / (n-k);
//     //Rcpp::Rcout << "ending" << std::endl;
//     //sigma = M2 / n - mu * mu.t();
//   } while ( norm(beta-beta_prev, 2) > eps && current_iter < max_iter);
//
//   // Last iteration
//   inv_sigma_lm = arma::inv_sympd(sigma_lm);
//   for(int i = 0; i < Y.n_rows; i++){
//     arma::mat mu_i =  X.row(i) * beta;
//     arma::mat N_posterior = c_posterior_approximation(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
//     arma::mat moments = c_moments_lrnm_hermite_precision(Y.row(i).t(),
//                                                          N_posterior.col(d), N_posterior.head_cols(d),
//                                                          mu_i.t(), sigma_lm,
//                                                          Binv, order);
//
//     H.row(i) = moments.col(d).t();
//   }
//   beta = arma::inv(X.t() * X) * X.t() * H;
//   //Rcpp::Rcout << beta;
//
//   return Rcpp::List::create(beta, sigma_lm, current_iter);
// }

//' @export
// [[Rcpp::export]]
Rcpp::List c_fit_lm_lrnm_hermite_centered(arma::mat Y, arma::mat B, arma::mat X, int order,
                                          double eps, int max_iter){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;
  arma::vec alpha = c_dm_fit_alpha(Y);
  arma::mat Binv = pinv(B).t();
  arma::mat P = arma::mat(Y);
  P.each_row() += alpha.t();

  arma::mat H = arma::log(P) * B;

  arma::mat beta = arma::inv(X.t() * X) * X.t() * H, beta_prev;
  arma::mat R = H - X * beta;
  arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;


  int current_iter = 0;
  do{
    current_iter++;
    //Rcpp::Rcout << "sigma:" << std::endl << sigma_lm;
    inv_sigma_lm = arma::inv_sympd(sigma_lm);
    //Rcpp::Rcout << "inv sigma:" << std::endl << inv_sigma_lm;
    beta_prev = arma::mat(beta);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);
    for(int i = 0; i < Y.n_rows; i++){
      arma::mat mu_i =  X.row(i) * beta;
      //Rcpp::Rcout << "start optimisation";
      //Rcpp::Rcout << mu_i << inv_sigma_lm;
      arma::mat N_posterior = c_posterior_approximation(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
      arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
                                                              N_posterior.col(d), N_posterior.head_cols(d),
                                                              mu_i.t(), sigma_lm,
                                                              Binv, order);
      H.row(i) = moments.col(d).t();
      M1 += moments.col(d);
      M2 += moments.head_cols(d);
    }
    //Rcpp::Rcout << "ending" << std::endl;
    beta = arma::inv(X.t() * X) * X.t() * H;
    //R = H - X * beta;
    sigma_lm = M2 / n; //R.t() * R / (n-k);
    //Rcpp::Rcout << "ending" << std::endl;
    //sigma = M2 / n - mu * mu.t();
  } while ( norm(beta-beta_prev, 2) > eps && current_iter < max_iter);

  // Last iteration
  inv_sigma_lm = arma::inv_sympd(sigma_lm);
  for(int i = 0; i < Y.n_rows; i++){
    arma::mat mu_i =  X.row(i) * beta;
    arma::mat N_posterior = c_posterior_approximation(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
    arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
                                                            N_posterior.col(d), N_posterior.head_cols(d),
                                                            mu_i.t(), sigma_lm,
                                                            Binv, order);

    H.row(i) = moments.col(d).t();
  }
  beta = arma::inv(X.t() * X) * X.t() * H;
  //Rcpp::Rcout << beta;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}

// /*
//  * Maybe following functions should be removed
//  */
//
// // [[Rcpp::export]]
// double c_d_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
//                          arma::mat Binv, int order = 100, int step_by = 100,
//                          double eps = 0.000001, int max_steps = 10){
//   double vcurrent = -1;
//   double vnext = c_d_lrnm_hermite(x, mu, sigma, Binv, order);
//   int step = 1;
//   while(abs(vcurrent - vnext) > eps &  step < max_steps){
//     step++;
//     order+=step_by;
//     vcurrent = vnext;
//     vnext = c_d_lrnm_hermite(x, mu, sigma, Binv, order);
//   }
//   return(vnext);
// }
//
// // [[Rcpp::export]]
// arma::vec c_m1_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
//                              arma::mat Binv, int order = 100, int step_by = 100,
//                              double eps = 0.000001, int max_steps = 10){
//   arma::vec vcurrent(mu.n_elem);
//   arma::vec vnext = c_m1_lrnm_hermite(x, mu, sigma, Binv, order);
//   int step = 1;
//   while(max(abs(vcurrent - vnext)) > eps &  step < max_steps){
//     step++;
//     order+=step_by;
//     vcurrent = vnext;
//     vnext = c_m1_lrnm_hermite(x, mu, sigma, Binv, order);
//   }
//   return(vnext);
// }
//
// // [[Rcpp::export]]
// arma::mat c_m2_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
//                              arma::mat Binv, int order = 100, int step_by = 100,
//                              double eps = 0.000001, int max_steps = 10){
//   arma::mat vcurrent(mu.n_elem, mu.n_elem);
//   arma::mat vnext = c_m2_lrnm_hermite(x, mu, sigma, Binv, order);
//   int step = 1;
//   while(max(max(abs(vcurrent - vnext))) > eps &  step < max_steps){
//     step++;
//     order+=step_by;
//     vcurrent = vnext;
//     vnext = c_m2_lrnm_hermite(x, mu, sigma, Binv, order);
//   }
//   return(vnext);
// }
//
// // [[Rcpp::export]]
// Rcpp::List c_fit_lrnm_hermite_(arma::mat X, arma::vec mu0, arma::mat sigma0,
//                                arma::mat Binv, int order = 100, int step_by = 100,
//                                double eps = 0.000001, int max_steps = 10,
//                                int em_max_steps = 10){
//   X = X.t();
//   int n = X.n_cols;
//
//   arma::mat H(mu0.n_elem, n);
//   arma::vec M1(mu0.n_elem);
//   arma::mat M2(mu0.n_elem, mu0.n_elem);
//   arma::vec mu = mu0;
//   arma::mat sigma = sigma0;
//   arma::vec mu_prev = mu0 + 1;
//
//
//   int step = 0;
//   while(max(abs(mu_prev - mu)) > eps & step < em_max_steps){
//     step++;
//     M1.zeros();
//     M2.zeros();
//     for(int i = 0; i < n; i++){
//       double prob = c_d_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps);
//       H.col(i) = c_m1_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps)/prob;
//       M1 += H.col(i);
//       M2 += (c_m2_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps)/prob);
//     }
//     mu_prev = mu;
//     mu = M1 / n;
//     sigma = M2/n - mu * mu.t();
//   }
//
//   return Rcpp::List::create(mu, sigma, H.t());
// }
