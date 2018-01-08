// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
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

// [[Rcpp::export]]
double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
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
arma::mat alr_basis(int K){
  arma::mat B = arma::zeros(K, K-1);
  B.row(K-1).fill(-1);
  B.diag().fill(1);
  //for(int i = 0; i< K-1; i++) B(i,i) = 1;
  return(B);
}

//' @export
// [[Rcpp::export]]
arma::vec mvf_maximum(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat B, double eps, int max_iter, double prop) {
  arma::mat A = alr_basis(B.n_rows).t();
  arma::mat Ainv = pinv(A);
  //Rcpp::Rcout << A;
  //Rcpp::Rcout << B;
  arma::vec mu_alr = A * B * mu_ilr;

  arma::mat sigma_alr = A * B * sigma_ilr * B.t() * A.t();
  arma::mat inv_sigma_alr = inv_sympd(sigma_alr);
  //Rcpp::Rcout << sigma_alr;

  int k = x.size() - 1;
  arma::vec a = arma::vec(k);

  for(int j = 0; j < k; j++){
    if(x[j] != 0 && x[k] != 0){
      a[j] = log(x[j]/x[k]);
    }else{
      if(x[j] == 0 && x[k] == 0){
        a[j] = 0;
      }else{
        if(x[j] == 0){
          a[j] = log(prop/x[k]);
        }else{
          a[j] = log(x[j]/prop);
        }
      }
    }
  }

  arma::vec out = arma::vec(a);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{
    current_iter++;
    for(int I=0; I<k; I++){
      deriv[I] =  mvf_deriv(I, out, mu_alr, inv_sigma_alr, x);
      for(int J=0; J<k; J++){
        deriv2(I,J) = mvf_deriv2(I, J, out, mu_alr, inv_sigma_alr, x);
      }
    }
    step = arma::solve(deriv2, deriv);
    out = out - step;
  }while( norm(step, 2) > eps && current_iter < max_iter);

  return B.t() * Ainv * out;
}
