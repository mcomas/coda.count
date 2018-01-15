// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
# include "nm.h"
# include "nm_expect.h"
# include "nm_optimise.h"
# include "nm_expect_mc_variations.h"
# include "coda.h"
#include "hermite.h"
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma){
  int k = H.n_cols;
  int n = H.n_rows;
  double log_det_val;
  double sign;
  log_det(log_det_val, sign, inv_sigma);
  double norm_const = -0.5 * k * log(2*PI) + 0.5 * log_det_val;
  arma::vec log_norm = arma::vec(n);

  for(int i = 0; i < n; i++){
    arma::mat cal = (H.row(i) - mu.t()) * inv_sigma * (H.row(i) - mu.t()).t();
    log_norm[i] = -0.5 * cal(0);
  }
  arma::vec norm = norm_const + log_norm;

  return(norm);
}

//' @export
// [[Rcpp::export]]
double lpmultinomial_const(arma::vec x){
  int K = x.size();
  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];
  arma::vec log_sums = arma::vec(1+x_total);
  log_sums[0] = 0;
  for(int l = 1; l <= x_total; l++) log_sums[l] = log_sums[l-1] + log(l);
  double constant = log_sums[x_total];
  for(int i = 0; i < K; i++) constant -= log_sums[x[i]];
  return(constant);
}


double lpmultinomial(arma::vec x, arma::vec p, double lconst){
  //double lconst = lpmultinomial_const(x);
  return( lconst + arma::dot(log(p),x) );
}

arma::mat lpnm_join(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::mat P, arma::mat H){
  double lconst = lpmultinomial_const(x);
  arma::mat lmult = arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1);
  arma::mat lnormal = ldnormal(H, mu, inv_sigma);
  return(lconst + lmult + lnormal);
}

arma::vec lpmultinomial_mult(arma::mat P, arma::vec x){
  return(arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1) );
}

double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                             arma::vec p, arma::vec h){
  double lmult = arma::accu(x % log(p));
  arma::vec y = h-mu;
  double lnormal = -0.5 * ((arma::mat)(y.t() * inv_sigma * y))(0,0);
  return(lmult + lnormal);
}

// [[Rcpp::export]]
double c_dlrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  double vcurrent = -1;
  double vnext = c_dnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(abs(vcurrent - vnext) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_dnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}


arma::vec c_m1_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  arma::vec vcurrent(mu.n_elem);
  arma::vec vnext = m1_lrnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(max(abs(vcurrent - vnext)) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = m1_lrnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}


arma::mat c_m2_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  arma::mat vcurrent(mu.n_elem, mu.n_elem);
  arma::mat vnext = m2_lrnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(max(max(abs(vcurrent - vnext))) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = m2_lrnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_hermite(arma::mat X, arma::vec mu0, arma::mat sigma0,
                   int order = 100, int step_by = 100,
                   double eps = 0.000001, int max_steps = 10,
                   int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;

  arma::mat H(mu0.n_elem, n);
  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;
  arma::vec mu_prev = mu0 + 1;


  int step = 0;
  while(max(abs(mu_prev - mu)) > eps & step < em_max_steps){
    step++;
    M1.zeros();
    M2.zeros();
    for(int i = 0; i < n; i++){
      double prob = c_dlrnm_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps);
      H.col(i) = c_m1_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps)/prob;
      M1 += H.col(i);
      M2 += (c_m2_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps)/prob);
    }
    mu_prev = mu;
    mu = M1 / n;
    sigma = M2/n - mu * mu.t();
  }

  return Rcpp::List::create(mu, sigma, inv_ilr_coordinates(H.t()));
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_maximum(arma::mat X, arma::vec mu0, arma::mat sigma0,
                              double tol = 10e-6, int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;
  int K = X.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat H(mu0.n_elem, n);
  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;
  arma::vec mu_prev = mu0 + 1;

  int step = 0;
  while(max(abs(mu_prev - mu)) > tol & step < em_max_steps){
    // Rcpp::Rcout << mu << std::endl;
    step++;
    for(int i = 0; i < n; i++){
      H.col(i) = mvf_maximum(X.col(i), mu, sigma, B, 10e-6, 50, 0.66);
    }
    mu_prev = mu;
    mu = mean(H,1);
    sigma = cov(H.t());
  }
  return Rcpp::List::create(mu, sigma, inv_ilr_coordinates(H.t()));
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_mc_init(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z,
                              double tol = 0.001, int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;
  int K = X.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat H(mu0.n_elem, n);

  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);

  arma::mat M1i(mu0.n_elem, n);
  arma::cube M2i(mu0.n_elem, mu0.n_elem, n);

  arma::vec mu_prev = mu0 + 1;
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;

  arma::mat Hs(Z.n_rows, Z.n_cols);

  M1.zeros();
  M2.zeros();
  for(int i = 0; i < n; i++){
    arma::vec x = X.col(i);
    arma::vec sampling_mu = mvf_maximum(x, mu, sigma, B, 0.0001, 50, 0.66);
    arma::vec lik_std = expected_mc_01_init(x, mu, sigma, Z, sampling_mu, Hs);
    arma::vec mu_exp  = expected_mc_mean(x, Hs, lik_std);

    M1 += mu_exp;
    M2 += expected_mc_var(x, mu, Hs, lik_std);

    M1i.col(i) = mu_exp;
    M2i.slice(i) = expected_mc_var(x, mu_exp, Hs, lik_std);
  }
  mu_prev = mu;
  mu = M1 / n;
  sigma = M2 / n;

  return Rcpp::List::create(mu, sigma, M1i, M2i);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_mc_step(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z,
                              arma::mat M1i, arma::cube M2i,
                              double tol = 0.001, int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;
  int K = X.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat H(mu0.n_elem, n);

  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);

  arma::vec mu_prev = mu0 + 1;
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;

  arma::mat Hs(Z.n_rows, Z.n_cols);

  M1.zeros();
  M2.zeros();
  for(int i = 0; i < n; i++){
    arma::vec x = X.col(i);
    arma::vec sampling_mu = M1i.col(i);
    arma::mat sampling_sigma = M2i.slice(i);


    //if(i == 0){

    //}

    if(min(eig_sym(sampling_sigma)) > 10e-15){
      //Rcpp::Rcout << "x" << std::endl << x << std::endl;
      //Rcpp::Rcout << "highest " << max(eig_sym(sampling_sigma) ) << " lowest eigenvalue " << min(eig_sym(sampling_sigma) ) << std::endl;
      arma::vec lik_std = expected_mc_03_init(x, mu, sigma, Z, sampling_mu, sampling_sigma, Hs);
      /*
      Rcpp::Rcout << "mu" << std::endl << mu << std::endl;
      Rcpp::Rcout << "sigma" << std::endl << sigma << std::endl;
      Rcpp::Rcout << "x" << std::endl << x << std::endl;
      Rcpp::Rcout << "sampling_mu" << std::endl << sampling_mu << std::endl;
      Rcpp::Rcout << "sampling_sigma" << std::endl << sampling_sigma << std::endl;
      //Rcpp::Rcout << "trace: " << trace(sampling_sigma) << std::endl;
      //Rcpp::Rcout << "lowest eigenvalue " << min(eig_sym(sampling_sigma) ) << std::endl;
      Rcpp::Rcout << "Hs" << std::endl << Hs << std::endl;
      Rcpp::Rcout << "Hs (mean)" << std::endl << mean(Hs,0) << std::endl;
      Rcpp::Rcout << "Hs (cov)" << std::endl << cov(Hs) << std::endl;
      */

      arma::vec mu_exp  = expected_mc_mean(x, Hs, lik_std);
      //Rcpp::Rcout << "Calculated" << std::endl;
      M1 += mu_exp;
      M2 += expected_mc_var(x, mu, Hs, lik_std);

      M1i.col(i) = mu_exp;
      M2i.slice(i) = expected_mc_var(x, mu_exp, Hs, lik_std);
    }else{
      M1 += M1i.col(i);
      M2 += (M1i.col(i) - mu) * (M1i.col(i) - mu).t();
    }
  }
  mu_prev = mu;
  mu = M1 / n;
  sigma = M2 / n;

  return Rcpp::List::create(mu, sigma, M1i, M2i);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_mc(arma::mat X, arma::vec mu0, arma::mat sigma0, arma::mat Z,
                         double tol = 0.001, int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;
  int K = X.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat H(mu0.n_elem, n);

  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);

  arma::mat M1i(mu0.n_elem, n);
  arma::cube M2i(mu0.n_elem, mu0.n_elem, n);

  arma::vec mu_prev = mu0 + 1;
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;

  arma::mat Hs(Z.n_rows, Z.n_cols);

  M1.zeros();
  M2.zeros();
  for(int i = 0; i < n; i++){
    arma::vec x = X.col(i);
    arma::vec sampling_mu = mvf_maximum(x, mu, sigma, B, 0.0001, 50, 0.66);
    arma::vec lik_std = expected_mc_01_init(x, mu, sigma, Z, sampling_mu, Hs);
    arma::vec mu_exp  = expected_mc_mean(x, Hs, lik_std);

    M1 += mu_exp;
    M2 += expected_mc_var(x, mu, Hs, lik_std);

    M1i.col(i) = mu_exp;
    M2i.slice(i) = expected_mc_var(x, mu_exp, Hs, lik_std);
  }
  mu_prev = mu;
  mu = M1 / n;
  sigma = M2 / n;

  // Starting iterations
  int step = 0;
  while(max(abs(mu_prev - mu)) > tol & step < em_max_steps){
    //Rcpp::Rcout << "Mu prev: " << mu_prev.t() << std::endl;
    Rcpp::Rcout << "Mu " << mu.t() << std::endl;
    Rcpp::checkUserInterrupt();
    step++;
    Rcpp:Rcout << "Step " << step << ". Abs. dif" << max(abs(mu_prev - mu)) << std::endl;
    M1.zeros();
    M2.zeros();
    arma::mat inv_sigma = inv_sympd(sigma);
    for(int i = 0; i < n; i++){
      arma::vec x = X.col(i);
      arma::vec sampling_mu = M1i.col(i);
      arma::mat sampling_sigma = M2i.slice(i);

      //Rcpp::Rcout << "Total:" << accu(x) << ", max eigen: " << max(eig_sym(sampling_sigma)) << ", min eigen: " << min(eig_sym(sampling_sigma)) << std::endl;
      arma::vec eigenval = eig_sym(sampling_sigma);

      if(max(eigenval) > 1e-10){
        //Rcpp::Rcout << "x" << std::endl << x << std::endl;
        //Rcpp::Rcout << "highest " << max(eig_sym(sampling_sigma) ) << " lowest eigenvalue " << min(eig_sym(sampling_sigma) ) << std::endl;

        arma::vec lik_std = expected_mc_03_init(x, mu, inv_sigma, Z, sampling_mu, sampling_sigma, Hs);

        /*
        Rcpp::Rcout << "mu" << std::endl << mu << std::endl;
        Rcpp::Rcout << "sigma" << std::endl << sigma << std::endl;
        Rcpp::Rcout << "x" << std::endl << x << std::endl;
        Rcpp::Rcout << "sampling_mu" << std::endl << sampling_mu << std::endl;
        Rcpp::Rcout << "sampling_sigma" << std::endl << sampling_sigma << std::endl;
        //Rcpp::Rcout << "trace: " << trace(sampling_sigma) << std::endl;
        //Rcpp::Rcout << "lowest eigenvalue " << min(eig_sym(sampling_sigma) ) << std::endl;
        Rcpp::Rcout << "Hs" << std::endl << Hs << std::endl;
        Rcpp::Rcout << "Hs (mean)" << std::endl << mean(Hs,0) << std::endl;
        Rcpp::Rcout << "Hs (cov)" << std::endl << cov(Hs) << std::endl;
         */

        arma::vec mu_exp  = expected_mc_mean(x, Hs, lik_std);

        M1 += mu_exp;
        M2 += expected_mc_var(x, mu, Hs, lik_std);

        M1i.col(i) = mu_exp;
        M2i.slice(i) = expected_mc_var(x, mu_exp, Hs, lik_std);

      }else{
        M1 += M1i.col(i);
        M2 += (M1i.col(i) - mu) * (M1i.col(i) - mu).t();
      }
    }
    mu_prev = mu;
    mu = M1 / n;
    sigma = M2 / n;
  }

  return Rcpp::List::create(mu, sigma, inv_ilr_coordinates(M1i.t()));
}
