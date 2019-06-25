// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::vec ldnormal2(arma::mat H, arma::vec mu, arma::mat inv_sigma){
  int k = H.n_cols;
  int n = H.n_rows;
  double log_det_val;
  double sign;
  log_det(log_det_val, sign, inv_sigma);
  double norm_const = -0.5 * k * log2pi + 0.5 * log_det_val;
  arma::vec log_norm = arma::vec(n);
  arma::mat x;
  for(int i = 0; i < n; i++){
    x = H.row(i) - mu.t();
    x = x * inv_sigma * x.t();
    log_norm[i] = -0.5 * x(0);
  }
  arma::vec norm = norm_const + log_norm;

  return(norm);
}


//' @export
// [[Rcpp::export]]
arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma){
  int k = H.n_cols;
  int n = H.n_rows;
  double log_det_val;
  double sign;
  log_det(log_det_val, sign, inv_sigma);
  double norm_const = -0.5 * k * log2pi + 0.5 * log_det_val;
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

//' @export
// [[Rcpp::export]]
arma::mat lpnm_join(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::mat P, arma::mat H){
  double lconst = lpmultinomial_const(x);
  arma::mat lmult = arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1);
  arma::mat lnormal = ldnormal(H, mu, inv_sigma);
  return(lconst + lmult + lnormal);
}

//' @export
// [[Rcpp::export]]
double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                             arma::vec p, arma::vec h){
  double lmult = arma::accu(x % log(p));
  arma::vec y = h-mu;
  double lnormal = -0.5 * ((arma::mat)(y.t() * inv_sigma * y))(0,0);
  return(lmult + lnormal);
}

//' @export
// [[Rcpp::export]]
double lpnm_join_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  arma::mat log_norm =  -(a-mu).t() * inv_sigma(arma::span::all, I);
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = 0;
  mult += x(I) * (kappa-exp(a(I))) / kappa;
  for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
  for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
  return log_norm(0) + mult;
}

//' @export
// [[Rcpp::export]]
double lpnm_join_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = -inv_sigma(I, J);
  if(I == J){
    for(int i = 0; i < k + 1; i++) mult -= x(i) * exp(a(I)) * (kappa-exp(a(I))) / (kappa*kappa);
  }else{
    for(int i = 0; i < k + 1; i++) mult += x(i) * exp(a(I) + a(J)) / (kappa*kappa);
  }
  return mult;
}

//' @export
// [[Rcpp::export]]
arma::vec lpnm_join_maximum_alr(arma::vec x, arma::vec mu_alr, arma::mat& inv_sigma_alr,
                                arma::vec a, double eps, int max_iter) {

  int k = x.size() - 1;

  arma::vec out = arma::vec(a);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{
    current_iter++;
    for(int I=0; I<k; I++){
      deriv[I] =  lpnm_join_deriv(I, out, mu_alr, inv_sigma_alr, x);
      for(int J=0; J<k; J++){
        deriv2(I,J) = lpnm_join_deriv2(I, J, out, mu_alr, inv_sigma_alr, x);
      }
    }
    step = arma::solve(deriv2, deriv);
    out = out - step;
  }while( norm(step, 2) > eps && current_iter < max_iter);

  return out;
}
