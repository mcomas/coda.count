// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "lrnm_gaussian_approx.h"
#include "lrnm_utils.h"

// //' @export
// // [[Rcpp::export]]
// arma::mat c_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma,
//                                         arma::mat &Binv, double eps = 1e-05, int niter = 1000){
//   unsigned d = Binv.n_cols;
//
//   arma::mat N_posterior(d,d+1);
//   // Rcpp::Rcout << inv_sigma << std::endl;
//   // Rcpp::Rcout << mu << std::endl;
//   // Rcpp::Rcout << x << std::endl;
//   N_posterior.col(d) = l_lrnm_join_maximum(x, mu, inv_sigma, Binv, eps, niter);
//   arma::mat D2 = l_lrnm_join_d2(N_posterior.col(d), x, mu, inv_sigma, Binv);
//   // Rcpp::Rcout << N_posterior.col(d);
//   // Rcpp::Rcout << -D2;
//   N_posterior.head_cols(d) = arma::inv_sympd(-D2);
//   return(N_posterior);
// }

//' @export
// [[Rcpp::export]]
arma::mat c_posterior_approximation_vec_sigma_inverse(arma::vec x, arma::vec mu,
                                                      arma::mat &inv_sigma,
                                                      arma::mat &Binv){
  unsigned d = Binv.n_cols;

  arma::mat N_posterior(d,d+1);
  N_posterior.col(d) = l_lrnm_join_maximum(x, mu, inv_sigma, Binv);
  arma::mat D2 = l_lrnm_join_d2(N_posterior.col(d), x, mu, inv_sigma, Binv);

  N_posterior.head_cols(d) = -D2;
  if(!N_posterior.head_cols(d).is_sympd()){
    Rcpp::Rcout << N_posterior.col(d) << std::endl <<
      x << std::endl << mu << std::endl << inv_sigma;
  }
  return(N_posterior);
}

// //' @export
// // [[Rcpp::export]]
// arma::mat c_lrnm_cond_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2,
//                                                   arma::mat &Binv, double eps = 1e-05, int niter = 1000){
//   unsigned d = Binv.n_cols - h2.size();
//
//   arma::mat N_posterior(d,d+1);
//   N_posterior.col(d) = l_lrnm_cond_join_maximum(x, mu, inv_sigma, h2, Binv, eps, niter);
//   arma::mat D2 = l_lrnm_cond_join_d2(N_posterior.col(d), x, mu, inv_sigma, h2, Binv);
//   N_posterior.head_cols(d) = arma::inv_sympd(-D2, arma::inv_opts::allow_approx);
//   return(N_posterior);
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_lrnm_cond_posterior_approximation_vec_sigma_inverse(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2,
//                                                                 arma::mat &Binv, double eps = 1e-05, int niter = 1000){
//   unsigned d = Binv.n_cols - h2.size();
//
//   arma::mat N_posterior(d,d+1);
//   // Rcpp::Rcout << inv_sigma << std::endl;
//   // Rcpp::Rcout << mu << std::endl;
//   // Rcpp::Rcout << x << std::endl;
//   N_posterior.col(d) = l_lrnm_cond_join_maximum(x, mu, inv_sigma, h2, Binv, eps, niter);
//   arma::mat D2 = l_lrnm_cond_join_d2(N_posterior.col(d), x, mu, inv_sigma, h2, Binv);
//   // Rcpp::Rcout << N_posterior.col(d);
//   // Rcpp::Rcout << -D2;
//   N_posterior.head_cols(d) = -D2;
//   return(N_posterior);
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_indep_lrnm_cond_posterior_approximation_vec(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2,
//                                                         arma::mat &Binv, double eps = 1e-05, int niter = 1000){
//
//   arma::mat N_posterior(1,2);
//   // Rcpp::Rcout << inv_sigma << std::endl;
//   // Rcpp::Rcout << mu << std::endl;
//   // Rcpp::Rcout << x << std::endl;
//   N_posterior.col(1) = l_indep_lrnm_cond_join_maximum(h1, ih1, x, mu, inv_sigma, h2, Binv, eps, niter);
//   arma::vec h(h1);
//   h(ih1) = N_posterior(0,1);
//   Rcpp::Rcout << h << std::endl;
//   arma::mat D2 = l_indep_lrnm_cond_join_d2(h, ih1, x, mu, inv_sigma, h2, Binv);
//   // Rcpp::Rcout << N_posterior.col(d);
//   // Rcpp::Rcout << -D2;
//   N_posterior.col(0) = -1/D2(0,0);
//   return(N_posterior);
// }
//
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_posterior_approximation_vec2(arma::vec x, arma::vec mu, arma::mat &inv_sigma,
// //                                          arma::mat &Binv, double eps = 1e-05, int niter = 1000){
// //   unsigned d = Binv.n_cols;
// //
// //   /*
// //    * https://stackoverflow.com/questions/48348079/applying-the-optim-function-in-r-in-c-with-rcpp
// //    * https://stackoverflow.com/questions/46473895/calling-rs-optim-function-from-within-c-using-rcpp
// //    *
// //    */
// //   arma::mat N_posterior(d,d+1);
// //
// //   // Extract R's optim function
// //   Rcpp::Environment stats("package:stats");
// //   Rcpp::Function optim = stats["optim"];
// //
// //   // Call the optim function from R in C++
// //   Rcpp::List opt_results = optim(Rcpp::_["par"]    = mu,
// //                                  // Make sure this function is not exported!
// //                                  Rcpp::_["fn"]     = Rcpp::InternalFunction(&neg_l_lrnm_join_no_constant_vec),
// //                                  Rcpp::_["method"] = "BFGS",
// //                                  // Pass in the other parameters as everything
// //                                  // is scoped environmentally
// //                                  Rcpp::_["x"] = x,
// //                                  Rcpp::_["mu"] = mu,
// //                                  Rcpp::_["inv_sigma"] = inv_sigma,
// //                                  Rcpp::_["Binv"] = Binv);
// //
// //   // Extract out the estimated parameter values
// //   arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
// //
// //   // Return estimated values
// //   return out;
// // }
//
// // //' @export
// // // [[Rcpp::export]]
// // arma::mat c_posterior_approximation_vec3(arma::vec x, arma::vec mu, arma::mat &inv_sigma,
// //                                          arma::mat &Binv, double eps = 1e-05, int niter = 1000){
// //   unsigned d = Binv.n_cols;
// //
// //   /*
// //    * https://stackoverflow.com/questions/48348079/applying-the-optim-function-in-r-in-c-with-rcpp
// //    * https://stackoverflow.com/questions/46473895/calling-rs-optim-function-from-within-c-using-rcpp
// //    *
// //    */
// //   arma::mat N_posterior(d,d+1);
// //
// //   // Extract R's optim function
// //   Rcpp::Environment stats("package:stats");
// //   Rcpp::Function optim = stats["optim"];
// //
// //   // Call the optim function from R in C++
// //   Rcpp::List opt_results = optim(Rcpp::_["par"]    = mu,
// //                                  // Make sure this function is not exported!
// //                                  Rcpp::_["fn"]     = Rcpp::InternalFunction(&neg_l_lrnm_join_no_constant_vec),
// //                                  Rcpp::_["gr"]     = Rcpp::InternalFunction(&neg_l_lrnm_join_d1),
// //                                  Rcpp::_["method"] = "BFGS",
// //                                  // Pass in the other parameters as everything
// //                                  // is scoped environmentally
// //                                  Rcpp::_["x"] = x,
// //                                  Rcpp::_["mu"] = mu,
// //                                  Rcpp::_["inv_sigma"] = inv_sigma,
// //                                  Rcpp::_["Binv"] = Binv);
// //
// //   // Extract out the estimated parameter values
// //   arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
// //
// //   // Return estimated values
// //   return out;
// // }
//
// //' @export
// // [[Rcpp::export]]
// arma::cube c_posterior_approximation(arma::mat X, arma::vec mu, arma::mat &sigma, arma::mat &B,
//                                      double eps = 1e-05, int niter = 1000){
//   arma::mat Binv = pinv(B).t();
//   arma::mat inv_sigma = inv_sympd(sigma);
//   arma::mat Xt = X.t();
//   arma::cube approx = arma::cube(mu.n_elem, mu.n_elem+1, X.n_rows);
//   for(int i = 0; i < X.n_rows; i++){
//     approx.slice(i) = c_posterior_approximation_vec(Xt.col(i), mu, inv_sigma, Binv, eps, niter);
//   }
//   return(approx);
// }
//
// // [[Rcpp::export]]
// Rcpp::List c_fit_lrnm_gaussian_approx(arma::mat X, arma::mat B,
//                                       double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){
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
//     arma::mat inv_sigma = arma::inv_sympd(sigma);
//     current_iter++;
//
//     mu_prev = arma::vec(mu);
//     arma::vec M1 = arma::zeros(d);
//     arma::mat M2 = arma::zeros(d, d);
//     for(int i = 0; i < H.n_rows; i++){
//       arma::mat N_posterior = c_posterior_approximation_vec(X.row(i).t(), mu, inv_sigma, Binv);
//       M1 += N_posterior.col(d);
//       M2 += (N_posterior.head_cols(d) + N_posterior.col(d) * N_posterior.col(d).t());
//     }
//     mu = M1 / n;
//     sigma = M2 / n - mu * mu.t();
//   } while ( norm(mu-mu_prev, 2) > eps && current_iter < max_iter);
//
//   return Rcpp::List::create(mu, sigma, current_iter);
// }
