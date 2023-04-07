// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "dm.h"
#include "hermite.h"
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec c_d_lrnm_hermite_mat(arma::mat X,
                               arma::vec mu_prior, arma::mat sigma_prior,
                               arma::mat Binv, int order){

  unsigned d = X.n_rows - 1;
  unsigned n = X.n_cols;


  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
  arma::vec dens = arma::zeros(n);
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  for(int i=0;i<n;i++){
    arma::mat N_posterior = c_lrnm_posterior_approximation_vec(X.col(i), mu_prior, inv_sigma_prior, Binv, 1e-5, 1000);
    arma::vec mu = N_posterior.col(d);
    arma::mat sigma = N_posterior.head_cols(d);
    arma::mat inv_sigma = arma::inv_sympd(sigma);


    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, sigma);
    arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

    unsigned int index[d+1];
    for(unsigned int j = 0; j <= d; j++) index[j] = 0;
    int position = 0, k = 0;

    double l_cmult = l_multinomial_const(X.col(i));

    do{
      double w = 0;
      arma::vec h(d);
      for(unsigned int j = 0; j < d; j++){
        h(j) = uni_hermite(index[j],0);
        w += uni_hermite(index[j],1);
      }
      h = mu + rotation * h;
      arma::vec p = exp(Binv * h);
      dens(i) += exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
        l_dnormal_vec(h, mu, inv_sigma) +
        l_multinomial(X.col(i), p/accu(p), l_cmult));
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
  }
  return dens;
}

// [[Rcpp::export]]
List c_lrnm_posterior_moments_hermite(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma, int order){

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;
  int ds = d;

  arma::vec clr_eigval;
  arma::mat clr_eigvec;
  eig_sym(clr_eigval, clr_eigvec, clr_sigma);
  arma::mat B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
  ds = B.n_rows;

  arma::vec B_mu = B * clr_mu;
  arma::mat B_sigma = B * clr_sigma * B.t();
  arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::mat Bs_N_mu = arma::zeros(d,n);
  arma::cube Bs_N_sigma = arma::zeros(d,d,n);

  arma::mat clr_E1(D,n);
  arma::cube clr_E2(D,D,n);

  int current_iter = 0;
  arma::vec h(ds);
  arma::vec deriv1(ds);
  arma::mat deriv2(ds, ds);
  arma::vec step(ds);

  for(int k = 0; k<n; k++){
    double eps = 1e-5;
    int max_iter = 100;
    h = B_mu;
    do{
      current_iter++;

      arma::vec Bh = B.t() * h;
      Bh = Bh - Bh.max();

      arma::vec eBh = exp(Bh);

      double  w = arma::accu(eBh);

      arma::vec  wBi = B * eBh;
      deriv1 = -inv_B_sigma  * (h-B_mu) + B * X.col(k) - sum(X.col(k)) * wBi / w;

      arma::mat wBij(ds,ds);
      for(int i=0; i<ds; i++){
        for(int j=0; j<ds; j++){
          wBij(i,j) = arma::accu(B.row(i).t() % B.row(j).t() % eBh);
        }
      }
      deriv2 = -inv_B_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

      step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
      h = h - 0.9 * step;

    }while( norm(step, 2) > eps && current_iter < max_iter);
    arma::vec B_N_mu = h;
    arma::mat B_N_inv_sigma = -deriv2;
    arma::mat B_N_sigma = arma::inv_sympd(B_N_inv_sigma);

    Bs_N_mu.col(k).subvec(0, ds-1) = B_N_mu;
    Bs_N_sigma.slice(k).submat(0, 0, ds-1, ds-1) = B_N_sigma;

    double l_cmult = l_multinomial_const(X.col(k));


    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, B_N_sigma);
    arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

    unsigned int index[ds+1];
    for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
    int position = 0, kappa = 0;
    arma::vec l_dens(std::pow(order, ds));
    double l_dens_max = 0;
    do{
      double w = 0;
      for(unsigned int i = 0; i < ds; i++){
        h(i) = uni_hermite(index[i],0);
        w += uni_hermite(index[i],1);
      }
      h = B_N_mu + rotation * h;
      arma::vec Bh = B.t() * h;
      Bh = Bh - Bh.max();

      arma::vec p = exp(Bh);
      l_dens(kappa) = w + l_dnormal_vec(h, B_mu, inv_B_sigma) -
        l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
        arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

      if(l_dens(kappa) > l_dens_max) l_dens_max = l_dens(kappa);


      // Calculate next coordinate
      index[position]++;
      while(index[position] == order){
        index[position] = 0;
        position++;
        index[position]++;
      }
      position = 0;
      kappa++;
    } while (index[ds] == 0);
    arma::vec dens = exp(l_dens-l_dens_max);

    double M0 = 0;
    arma::vec M1 = arma::zeros(ds);
    arma::mat M2 = arma::zeros(ds,ds);

    for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
    position = 0;
    kappa = 0;
    do{
      for(unsigned int i = 0; i < ds; i++){
        h(i) = uni_hermite(index[i],0);
      }
      h = B_N_mu + rotation * h;

      M0 += dens(kappa);
      M1 += h * dens(kappa);
      M2 += h * h.t() * dens(kappa);

      // Calculate next coordinate
      index[position]++;
      while(index[position] == order){
        index[position] = 0;
        position++;
        index[position]++;
      }
      position = 0;
      kappa++;
    } while (index[ds] == 0);

    clr_E1.col(k) = B.t() * (M1/M0);
    clr_E2.slice(k) = B.t() * (M2/M0) * B;
  }
  List res;
  res["clr_E1"] = clr_E1;
  res["clr_E2"] = clr_E2;
  res["B"] = B;
  res["Bs_N_mu"] = Bs_N_mu;
  res["Bs_N_sigma"] = Bs_N_sigma;
  res["ds"] = ds;
  return(res);
}


// [[Rcpp::export]]
List c_lrnm_fit_hermite(arma::mat& X, int order, double em_eps, int em_max_iter){
  List dm_ini = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = dm_ini["alpha"];
  arma::mat P(X);
  P.each_col() += alpha;

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;
  int ds = d;

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  arma::mat B;
  arma::mat Bs_N_mu = arma::zeros(d,n);
  arma::cube Bs_N_sigma = arma::zeros(d,d,n);

  arma::mat clr_E1(D,n);
  arma::cube clr_E2(D,D,n);

  arma::vec clr_mu_new;
  arma::mat clr_sigma_new;

  int em_iter = 0;
  double em_current_eps;
  do{
    em_iter++;

    arma::vec clr_eigval;
    arma::mat clr_eigvec;
    eig_sym(clr_eigval, clr_eigvec, clr_sigma);
    B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
    ds = B.n_rows;

    arma::vec B_mu = B * clr_mu;
    arma::mat B_sigma = B * clr_sigma * B.t();
    arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);

    arma::vec h(ds);
    arma::vec deriv1(ds);
    arma::mat deriv2(ds, ds);
    arma::vec step(ds);

    for(int k = 0; k<n; k++){
      double eps = 1e-5;
      int max_iter = 100;
      h = B * clr_H.col(k);
      int laplace_iter = 0;
      do{
        laplace_iter++;

        arma::vec Bh = B.t() * h;
        Bh = Bh - Bh.max();

        arma::vec eBh = exp(Bh);

        double  w = arma::accu(eBh);

        arma::vec  wBi = B * eBh;
        deriv1 = -inv_B_sigma  * (h-B_mu) + B * X.col(k) - sum(X.col(k)) * wBi / w;

        arma::mat wBij(ds,ds);
        for(int i=0; i<ds; i++){
          for(int j=0; j<ds; j++){
            wBij(i,j) = arma::accu(B.row(i).t() % B.row(j).t() % eBh);
          }
        }
        deriv2 = -inv_B_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

        step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
        h = h - 0.9 * step;

      }while( norm(step, 2) > eps && laplace_iter < max_iter);
      arma::vec B_N_mu = h;
      arma::mat B_N_inv_sigma = -deriv2;
      arma::mat B_N_sigma = arma::inv_sympd(B_N_inv_sigma);

      Bs_N_mu.col(k).subvec(0, ds-1) = B_N_mu;
      Bs_N_sigma.slice(k).submat(0, 0, ds-1, ds-1) = B_N_sigma;


      double l_cmult = l_multinomial_const(X.col(k));

      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, B_N_sigma);
      arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

      unsigned int index[ds+1];
      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      int position = 0, kappa = 0;
      arma::vec l_dens(std::pow(order, ds));
      double l_dens_max = 0;
      do{
        double w = 0;
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
          w += uni_hermite(index[i],1);
        }
        h = B_N_mu + rotation * h;
        arma::vec Bh = B.t() * h;
        Bh = Bh - Bh.max();

        arma::vec p = exp(Bh);
        l_dens(kappa) = w + l_dnormal_vec(h, B_mu, inv_B_sigma) -
          l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
          arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

        if(l_dens(kappa) > l_dens_max) l_dens_max = l_dens(kappa);


        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);
      arma::vec dens = exp(l_dens-l_dens_max);

      double M0 = 0;
      arma::vec M1 = arma::zeros(ds);
      arma::mat M2 = arma::zeros(ds,ds);

      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      position = 0;
      kappa = 0;
      do{
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
        }
        h = B_N_mu + rotation * h;

        M0 += dens(kappa);
        M1 += h * dens(kappa);
        M2 += h * h.t() * dens(kappa);

        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);

      clr_E1.col(k) = B.t() * (M1/M0);
      clr_E2.slice(k) = B.t() * (M2/M0) * B;


      clr_E1.col(k) = B.t() * (M1/M0);
      clr_E2.slice(k) = B.t() * (M2/M0) * B;
    }
    clr_mu_new = mean(clr_E1, 1);
    em_current_eps = norm(clr_mu-clr_mu_new, 1);

    clr_mu = clr_mu_new;
    arma::mat M = mean(clr_E2, 2);
    clr_sigma = M - clr_mu_new * clr_mu_new.t();

  } while (em_iter < em_max_iter & em_current_eps > em_eps);

  List res;
  res["clr_mu"] = clr_mu;
  res["clr_sigma"] = clr_sigma;
  res["ds"] = ds;
  res["clr_E1"] = clr_E1;
  res["em_iter"] = em_iter;
  res["Bs"] = B;
  res["Bs_N_mu"] = Bs_N_mu;
  res["Bs_N_sigma"] = Bs_N_sigma;
  return(res);
}

// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_hermite(arma::vec x,
//                                  arma::vec mu, arma::mat sigma,
//                                  arma::vec mu_prior, arma::mat sigma_prior,
//                                  arma::mat Binv, int order,
//                                  arma::vec mu_centering){
//   unsigned d = Binv.n_cols;
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
//
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_moments_lrnm_hermite_sigma_inverse(arma::vec x,
// //                                                arma::vec mu, arma::mat inv_sigma,
// //                                                arma::vec mu_prior, arma::mat inv_sigma_prior,
// //                                                arma::mat Binv, int order,
// //                                                arma::vec mu_centering){
// //   unsigned d = Binv.n_cols;
// //   arma::mat uni_hermite = hermite(order);
// //   uni_hermite.col(1) = log(uni_hermite.col(1));
// //
// //   arma::vec eigval;
// //   arma::mat eigvec;
// //
// //   eig_sym(eigval, eigvec, inv_sigma);
// //   arma::mat rotation = eigvec * arma::diagmat(sqrt(1/eigval));
// //
// //   unsigned int index[d+1];
// //   for(unsigned int i = 0; i <= d; i++) index[i] = 0;
// //   int position = 0, k = 0;
// //   double M0 = 0;
// //   arma::vec M1 = arma::zeros(d);
// //   arma::mat M2 = arma::zeros(d,d);
// //   double l_cmult = l_multinomial_const(x);
// //
// //   do{
// //     double w = 0;
// //     arma::vec h(d);
// //     for(unsigned int i = 0; i < d; i++){
// //       h(i) = uni_hermite(index[i],0);
// //       w += uni_hermite(index[i],1);
// //     }
// //     h = mu + rotation * h;
// //     arma::vec p = exp(Binv * h);
// //     double dens = exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
// //                       l_dnormal_vec(h, mu, inv_sigma) +
// //                       l_multinomial(x, p/accu(p), l_cmult));
// //     M0 += dens;
// //     M1 += h * dens;
// //     M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
// //
// //     // Calculate next coordinate
// //     index[position]++;
// //     while(index[position] == order){
// //       index[position] = 0;
// //       position++;
// //       index[position]++;
// //     }
// //     position = 0;
// //     k++;
// //   } while (index[d] == 0);
// //
// //   arma::mat moments(d, d+1);
// //   moments.col(d) = M1/M0;
// //   moments.head_cols(d) = M2/M0;
// //   return moments;
// // }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_cond_hermite_1d(arma::vec x,
//                                          arma::vec mu, double sigma,
//                                          arma::vec Mc, arma::mat inv_Sc, arma::vec h2,
//                                          arma::mat Binv, int order,
//                                          double mu_centering){
//
//
//   arma::span I1 = arma::span(0, 0);
//   arma::span I2 = arma::span(1,  Binv.n_cols-1);
//
//   // arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);
//   //
//   // arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
//   // arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   // arma::mat inv_Sc = arma::inv_sympd(Sc, arma::inv_opts::allow_approx);
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
//   int position = 0, k = 0;
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(1);
//   arma::mat M2 = arma::zeros(1, 1);
//   double l_cmult = l_multinomial_const(x);
//   arma::vec h_mu = arma::vec(1+h2.size());
//   h_mu(I1) = mu;
//   h_mu(I2) = h2;
//   arma::vec p_mu = exp(Binv * h_mu);
//
//   double w;
//   arma::vec h1(1);
//   arma::vec h = arma::vec(1+h2.size());
//   h(I2) = h2;
//   for(int i = 0; i < order; i++){
//
//     h1 = mu + sqrt(sigma) * uni_hermite(i,0);
//     w  = uni_hermite(i,1);
//
//     h(I1) = h1;
//     arma::vec p = exp(Binv * h);
//     // if(false){
//     //   Rcpp::Rcout << "l_dnormal_vec(h1, Mc, inv_Sc):" << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
//     //   Rcpp::Rcout << "l_dnormal_vec(h1, mu, inv_sigma):" << l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
//     //   Rcpp::Rcout << "l_multinomial(x, p/accu(p), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
//     //   Rcpp::Rcout << "l_multinomial(x, p_mu/accu(p_mu), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
//     //   Rcpp::Rcout << "ldens: " << w + l_dnormal_vec(h1, Mc, inv_Sc) -
//     //     l_dnormal_vec(h1, mu, inv_sigma) +
//     //     l_multinomial(x, p/accu(p), l_cmult) -
//     //     l_multinomial(x, p_mu/accu(p_mu), l_cmult) << std::endl << std::endl;
//     // }
//     double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
//                       l_dnormal_vec(h1, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult) -
//                       l_multinomial(x, p_mu/accu(p_mu), l_cmult));  // For big x, low variability.
//
//
//     M0 += dens;
//     M1 += h1 * dens;
//     M2 += (h1-mu_centering) * (h1-mu_centering) * dens;
//     // Calculate next coordinate
//   }
//
//
//   arma::mat moments(1, 2);
//   moments.col(1) = M1/M0;
//   moments.head_cols(1) = M2/M0;
//   return moments;
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_moments_lrnm_cond_hermite(arma::vec x,
//                                       arma::vec mu, arma::mat sigma,
//                                       arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
//                                       arma::mat Binv, int order,
//                                       arma::vec mu_centering){
//   unsigned dh1 =  Binv.n_cols-h2.size();
//
//   arma::span I1 = arma::span(0, dh1-1);
//   arma::span I2 = arma::span(dh1,  Binv.n_cols-1);
//
//   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);
//
//   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
//   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
//   arma::mat inv_Sc = arma::inv_sympd(Sc, arma::inv_opts::allow_approx);
//
//   arma::mat uni_hermite = hermite(order);
//   uni_hermite.col(1) = log(uni_hermite.col(1));
//
//   arma::vec eigval;
//   arma::mat eigvec;
//
//   arma::mat inv_sigma = arma::inv_sympd(sigma, arma::inv_opts::allow_approx);
//
//   eig_sym(eigval, eigvec, sigma);
//   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
//
//   Rcpp::Rcout << rotation << std::endl;
//   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2), arma::inv_opts::allow_approx);
//   double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);
//
//   unsigned int index[dh1+1];
//   for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
//   int position = 0, k = 0;
//   double M0 = 0;
//   arma::vec M1 = arma::zeros(dh1);
//   arma::mat M2 = arma::zeros(dh1,dh1);
//   double l_cmult = l_multinomial_const(x);
//   arma::vec h_mu = join_cols(mu, h2);
//   arma::vec p_mu = exp(Binv * h_mu);
//   do{
//     double w = 0;
//     arma::vec h1(dh1);
//     for(unsigned int i = 0; i < dh1; i++){
//       h1(i) = uni_hermite(index[i],0);
//       w += uni_hermite(index[i],1);
//     }
//
//     h1 = mu + rotation * h1; // rotation 0x0
//
//     arma::vec h = join_cols(h1, h2);
//     arma::vec p = exp(Binv * h);
//     // Rcpp::Rcout << "l_ph2:" << l_ph2 << std::endl;
//     // Rcpp::Rcout << "l_dnormal_vec(h1, Mc, inv_Sc):" << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
//     // Rcpp::Rcout << "l_dnormal_vec(h1, mu, inv_sigma):" << l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
//     // Rcpp::Rcout << "l_multinomial(x, p/accu(p), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
//     // Rcpp::Rcout << "l_multinomial(x, p_mu/accu(p_mu), l_cmult)):" << l_multinomial(x, p/accu(p), l_cmult) << std::endl;
//     double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
//                       l_dnormal_vec(h1, mu, inv_sigma) +
//                       l_multinomial(x, p/accu(p), l_cmult) -
//                       l_multinomial(x, p_mu/accu(p_mu), l_cmult));  // For big x, low variability.
//
//     // Rcpp::Rcout << "ldens: " << w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
//     //   l_dnormal_vec(h1, mu, inv_sigma) +
//     //   l_multinomial(x, p/accu(p), l_cmult) << std::endl << std::endl;
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
//
//   arma::mat moments(dh1, dh1+1);
//   moments.col(dh1) = M1/M0;
//   moments.head_cols(dh1) = M2/M0;
//   return moments;
// }
// //
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_moments_lrnm_cond_hermite2(arma::vec x,
// //                                        arma::vec mu, arma::mat sigma,
// //                                        arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
// //                                        arma::mat Binv, int order,
// //                                        arma::vec mu_centering){
// //   unsigned dh1 =  Binv.n_cols-h2.size();
// //
// //   arma::span I1 = arma::span(0, dh1-1);
// //   arma::span I2 = arma::span(dh1,  Binv.n_cols-1);
// //
// //   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
// //
// //   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
// //   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
// //   arma::mat inv_Sc = arma::inv_sympd(Sc);
// //
// //   arma::mat uni_hermite = hermite(order);
// //   uni_hermite.col(1) = log(uni_hermite.col(1));
// //
// //   arma::vec eigval;
// //   arma::mat eigvec;
// //
// //   arma::mat inv_sigma = arma::inv_sympd(sigma);
// //
// //   eig_sym(eigval, eigvec, sigma);
// //   arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
// //
// //   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
// //   // double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);
// //
// //   unsigned int index[dh1+1];
// //   for(unsigned int i = 0; i <= dh1; i++) index[i] = 0;
// //   int position = 0, k = 0;
// //   double M0 = 0;
// //   arma::vec M1 = arma::zeros(dh1);
// //   arma::mat M2 = arma::zeros(dh1,dh1);
// //   double l_cmult = l_multinomial_const(x);
// //
// //   do{
// //     double w = 0;
// //     arma::vec h1(dh1);
// //     for(unsigned int i = 0; i < dh1; i++){
// //       h1(i) = uni_hermite(index[i],0);
// //       w += uni_hermite(index[i],1);
// //     }
// //     h1 = mu + rotation * h1;
// //     arma::vec h = join_cols(h1, h2);
// //     arma::vec p = exp(Binv * h);
// //     double dens = exp(w + l_dnormal_vec(h1, Mc, inv_Sc) -
// //                       l_dnormal_vec(h1, mu, inv_sigma) +
// //                       l_multinomial(x, p/accu(p), l_cmult));
// //     // Rcpp::Rcout << dens;
// //     M0 += dens;
// //     M1 += h1 * dens;
// //     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
// //     // Calculate next coordinate
// //     index[position]++;
// //     while(index[position] == order){
// //       index[position] = 0;
// //       position++;
// //       index[position]++;
// //     }
// //     position = 0;
// //     k++;
// //   } while (index[dh1] == 0);
// //
// //   arma::mat moments(dh1, dh1+1);
// //   moments.col(dh1) = M1/M0;
// //   moments.head_cols(dh1) = M2/M0;
// //   return moments;
// // }
// //
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_moments_lrnm_cond1_hermite(arma::vec h, int i1,
// //                                        arma::vec x, arma::vec mu, double sigma,
// //                                        arma::vec mu_prior, arma::mat sigma_prior,
// //                                        arma::mat Binv, int order,
// //                                        arma::vec mu_centering){
// //
// //
// //   arma::uvec I1(1);
// //   I1(0) = i1 - 1;
// //   arma::uvec I2(Binv.n_cols-1);
// //   for(int i = 0, k = 0; i < Binv.n_cols; i++){
// //     if(i + 1 != i1){
// //       I2(k) = i;
// //       k++;
// //     }
// //   }
// //   // Rcpp::Rcout << I1 << std::endl;
// //   // Rcpp::Rcout << I2 << std::endl;
// //
// //   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
// //
// //   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h(I2)-mu_prior(I2));
// //   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
// //   arma::mat inv_Sc = arma::inv_sympd(Sc);
// //
// //
// //   arma::mat uni_hermite = hermite(order);
// //   uni_hermite.col(1) = log(uni_hermite.col(1));
// //
// //   arma::mat inv_sigma(1,1);
// //   inv_sigma(0,0) = 1/sigma;
// //
// //   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
// //   double l_ph2 = l_dnormal_vec(h(I2), mu_prior(I2), inv_sigma_prior_2);
// //
// //   double M0 = 0;
// //   arma::vec M1 = arma::zeros(1);
// //   arma::mat M2 = arma::zeros(1,1);
// //   double l_cmult = l_multinomial_const(x);
// //
// //   // Rcpp::Rcout << uni_hermite << std::endl;
// //   double w = 0;
// //   arma::vec h1;
// // //   arma::vec h = arma::zeros(Binv.n_cols);
// // //   h(I2) = h2;
// //   for(int i = 0; i < order; i++){
// //     h1 = sqrt(sigma) * uni_hermite(i,0) + mu;
// //     w = uni_hermite(i,1);
// //     h(I1) = h1;
// //     arma::vec p = exp(Binv * h);
// //
// //     double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
// //                       l_dnormal_vec(h1, mu, inv_sigma) +
// //                       l_multinomial(x, p/accu(p), l_cmult));
// //
// //     M0 += dens;
// //     M1 += h1 * dens;
// //     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
// //   }
// //
// //
// //
// //   arma::mat moments(1, 2);
// //   moments.col(1) = M1/M0;
// //   moments.col(0) = M2/M0;
// //   return moments;
// // }
// //
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_moments_indep_lrnm_cond_hermite(arma::vec x, int i1,
// //                                             arma::vec mu, double sigma,
// //                                             arma::vec mu_prior, arma::mat sigma_prior, arma::vec h2,
// //                                             arma::mat Binv, int order,
// //                                             arma::vec mu_centering){
// //   unsigned dh1 =  Binv.n_cols-h2.size();
// //
// //   arma::uvec I1(1);
// //   I1(0) = i1 - 1;
// //   arma::uvec I2(Binv.n_cols-1);
// //
// //   for(int i = 0, k = 0; i < Binv.n_cols; i++){
// //     if(i + 1 != i1){
// //       I2(k) = i;
// //       k++;
// //     }
// //   }
// //
// //
// //   arma::mat inv_sigma_prior_2 = arma::inv_sympd(sigma_prior(I2,I2));
// //
// //   arma::mat Mc = mu_prior(I1) + sigma_prior(I1,I2) * inv_sigma_prior_2 * (h2-mu_prior(I2));
// //   arma::mat Sc = sigma_prior(I1,I1) - sigma_prior(I1,I2) * inv_sigma_prior_2 * sigma_prior(I2,I1);
// //   arma::mat inv_Sc = arma::inv_sympd(Sc);
// //
// //   arma::mat uni_hermite = hermite(order);
// //   uni_hermite.col(1) = log(uni_hermite.col(1));
// //
// //   arma::vec eigval;
// //   arma::mat eigvec;
// //
// //   arma::mat inv_sigma(1,1);
// //   inv_sigma(0,0) = 1/sigma;
// //
// //   arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior(I2,I2));
// //   double l_ph2 = l_dnormal_vec(h2, mu_prior(I2), inv_sigma_prior);
// //
// //   double M0 = 0;
// //   arma::vec M1 = arma::zeros(dh1);
// //   arma::mat M2 = arma::zeros(dh1,dh1);
// //   double l_cmult = l_multinomial_const(x);
// //
// //   //Rcpp::Rcout << uni_hermite << std::endl;
// //   double w = 0;
// //   arma::vec h1;
// //   arma::vec h = arma::zeros(Binv.n_cols);
// //   h(I2) = h2;
// //   for(int i = 0; i < order; i++){
// //     h1 = sqrt(sigma) * uni_hermite(i,0) + mu;
// //     w = uni_hermite(i,1);
// //     h(I1) = h1;
// //     arma::vec p = exp(Binv * h);
// //
// //     double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
// //                       l_dnormal_vec(h1, mu, inv_sigma) +
// //                       l_multinomial(x, p/accu(p), l_cmult));
// //
// //     M0 += dens;
// //     M1 += h1 * dens;
// //     M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
// //   }
// //
// //   /*
// //    do{
// //    double w = 0;
// //    arma::vec h1(dh1);
// //    for(unsigned int i = 0; i < dh1; i++){
// //    h1(i) = uni_hermite(index[i],0);
// //    w += uni_hermite(index[i],1);
// //    }
// //    h1 = mu + h1;
// //    arma::vec h = join_cols(h1, h2);
// //    arma::vec p = exp(Binv * h);
// //    // Rcpp::Rcout << "Nh2: " << l_ph2 << std::endl;
// //    // Rcpp::Rcout << "Nc: " << l_dnormal_vec(h1, Mc, inv_Sc) << std::endl;
// //    // Rcpp::Rcout << "N: " << -l_dnormal_vec(h1, mu, inv_sigma) << std::endl;
// //    // Rcpp::Rcout << "Mult: " << l_multinomial(x, p/accu(p), l_cmult) << std::endl << std::endl;
// //    double dens = exp(w + l_ph2 + l_dnormal_vec(h1, Mc, inv_Sc) -
// //    l_dnormal_vec(h1, mu, inv_sigma) +
// //    l_multinomial(x, p/accu(p), l_cmult));
// //    // Rcpp::Rcout << dens;
// //    M0 += dens;
// //    M1 += h1 * dens;
// //    M2 += (h1-mu_centering) * (h1-mu_centering).t() * dens;
// //    // Calculate next coordinate
// //    index[position]++;
// //    while(index[position] == order){
// //    index[position] = 0;
// //    position++;
// //    index[position]++;
// //    }
// //    position = 0;
// //    k++;
// //    } while (index[dh1] == 0);
// //    */
// //
// //   arma::mat moments(dh1, dh1+1);
// //   moments.col(dh1) = M1/M0;
// //   moments.head_cols(dh1) = M2/M0;
// //   return moments;
// // }
// //
// // // //' @export
// // // // [[Rcpp::export]]
// // // Rcpp::List c_obtain_moments_lrnm_hermite(arma::mat Y,
// // //                                          arma::vec mu, arma::mat sigma,
// // //                                          arma::mat B, int order){
// // //   int n = Y.n_rows;
// // //   int d = Y.n_cols - 1;
// // //
// // //   arma::mat Binv = pinv(B).t();
// // //   arma::mat inv_sigma = arma::inv_sympd(sigma);
// // //
// // //   arma::mat M1 = arma::zeros(d, n);
// // //   arma::cube M2 = arma::zeros(d, d, n);
// // //
// // //   for(int i = 0; i < Y.n_rows; i++){
// // //     arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu, inv_sigma, Binv);
// // //     arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
// // //                                                             N_posterior.col(d), N_posterior.head_cols(d),
// // //                                                             mu, sigma,
// // //                                                             Binv, order);
// // //     M1.col(i) = moments.col(d);
// // //     M2.slice(i) = moments.head_cols(d);
// // //   }
// // //   return Rcpp::List::create(M1, M2);
// // // }
// //
// // //' @export
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
//       arma::mat N_posterior = c_lrnm_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
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
//     arma::mat N_posterior = c_lrnm_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
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
