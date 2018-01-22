#include <RcppArmadillo.h>

#include "Eigen/Dense"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/bfgssolver.h"
#include "cppoptlib/solver/newtondescentsolver.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

template <typename MatrixType>
inline typename MatrixType::Scalar logdet(const MatrixType& M, bool use_cholesky = false) {
  using namespace Eigen;
  using std::log;
  typedef typename MatrixType::Scalar Scalar;
  Scalar ld = 0;
  if (use_cholesky) {
    LLT<Eigen::Matrix<Scalar,Dynamic,Dynamic>> chol(M);
    auto& U = chol.matrixL();
    for (unsigned i = 0; i < M.rows(); ++i)
      ld += log(U(i,i));
    ld *= 2;
  } else {
    PartialPivLU<Eigen::Matrix<Scalar,Dynamic,Dynamic>> lu(M);
    auto& LU = lu.matrixLU();
    Scalar c = lu.permutationP().determinant(); // -1 or 1
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const auto& lii = LU(i,i);
      if (lii < Scalar(0)) c *= -1;
      ld += log(abs(lii));
    }
    ld += log(c);
  }
  return ld;
}

template<typename T>
class ConditionalMode : public cppoptlib::Problem<T> {
public:
  using typename cppoptlib::Problem<T>::TVector;
  using typename cppoptlib::Problem<T>::THessian;
  using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

protected:
  TVector x;
  const TVector mu;
  const MatrixType inv_sigma;
  const double log2pi = std::log(2.0 * M_PI);

public:
  ConditionalMode(const TVector &mu_,
                  const MatrixType &inv_sigma_) : x(mu_.size()+1), mu(mu_), inv_sigma(inv_sigma_) {}

  void update_x(const TVector x_){
    for(int i = 0; i < x_.size(); i++){
      x(i) = x_(i);
    }
  }

  T value(const TVector &a) {

    int k = a.size();
    int K = k + 1;

    double norm_const = -0.5 * k * log2pi - 0.5 * logdet(inv_sigma, true);

    MatrixType log_norm =  -0.5 * (a-mu).transpose() * inv_sigma * (a-mu);

    double norm = norm_const + log_norm(0,0);

    double x_total = 0;
    for(int i = 0; i < K; i++) x_total += x(i);

    double x_total_fact = log(x_total);
    for(int l = 1; l < x_total; l++) x_total_fact += log(l);

    double x_parts_fact = 0;
    for(int i = 0; i < K; i++){
      for(int l = 1; l <= x(i); l++){
        x_parts_fact += log(l);
      }
    }

    double mult_const = 0;
    mult_const += (x_total_fact - x_parts_fact);

    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a[i]);
    double multinom = -x[k] * log(kappa);
    for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

    double mult = mult_const + multinom;
    return(-norm - mult);
  }

  void gradient(const TVector &a, TVector &grad){
    int k = a.size();
    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a(i));

    for(int I=0; I < k; I++){
      Eigen::MatrixXd log_norm = -(a-mu).transpose() * inv_sigma.col(I);
      double mult = 0;
      mult += x(I) * (kappa-exp(a(I))) / kappa;
      for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
      for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
      grad(I) = -log_norm(0) - mult;
    }
  }

  void hessian(const TVector &a, THessian &hessian) {
    int k = a.size();
    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a(i));

    for(int I = 0; I < k; I++){
      for(int J = 0; J < k; J++){
        hessian(I,J) = inv_sigma(I, J);
        if(I == J){
          for(int i = 0; i < k + 1; i++) hessian(I,J) += x(i) * exp(a(I)) * (kappa-exp(a(I))) / (kappa*kappa);
        }else{
          for(int i = 0; i < k + 1; i++) hessian(I,J) -= x(i) * exp(a(I) + a(J)) / (kappa*kappa);
        }
      }
    }
  }
};


// arma::vec bfgsSolver(arma::vec x, arma::vec mu, arma::mat sigma) {
//   int k = mu.n_elem;
//   Eigen::VectorXd x_(k+1);
//   Eigen::VectorXd mu_(k);
//   Eigen::MatrixXd sigma_(k, k);
//   for(int i=0; i<x.n_elem; i++) x_(i) = x(i);
//   for(int i=0; i<mu.n_elem; i++) mu_(i) = mu(i);
//   for(int i=0; i<sigma.n_rows; i++) for(int j=0; j<sigma.n_cols; j++) sigma_(i, j) = sigma(i, j);
//
//
//   typedef ConditionalMode<double> TConditionalMode;
//   typedef typename TConditionalMode::TVector TVector;
//   typedef typename TConditionalMode::MatrixType MatrixType;
//
//   TConditionalMode f(mu_, sigma_);
//   cppoptlib::BfgsSolver<TConditionalMode> solver;
//   Eigen::VectorXd x0(k);
//   for(int i=0; i < k; i++){
//     x0(i) = log(x(i)+0.5) - log(x(k)+0.5);
//   }
//   for(int i = 0; i < n; i++){
//     f.update(x);
//   }
//   solver.minimize(f, x0);
//
//   arma::vec res = arma::vec(k);
//   for(int i=0; i<k; i++) res(i) = x0(i);
//
//   return res;
// }


//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_maximum_alr(arma::mat X, arma::vec mu0, arma::mat sigma0,
                                  double tol = 10e-6, int em_max_steps = 100){

  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  Eigen::MatrixXd X_(K, n);
  for(int i = 0; i < K; i++) for(int j = 0; j < n; j++) X_(i,j) = X(j,i);

  Eigen::VectorXd mu0_(k);
  for(int i = 0; i < k; i++) mu0_(i) = mu0(i);

  Eigen::MatrixXd sigma0_(k, k);
  for(int i = 0; i < k; i++) for(int j = 0; j < k; j++) sigma0_(i, j) = sigma0(i, j);

  Rcpp::Rcout << X_ << std::endl;

  Eigen::MatrixXd A_(k, n);
  for(int i = 0; i < n; i++){
    double denom = log(X_(k,i) + 0.5);
    for(int j = 0; j < k; j++){
      A_(j,i) = log(X_(j,i) + 0.5) - denom;
    }
  }
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << A_ << std::endl;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
  Eigen::VectorXd mu(mu0_);
  Eigen::MatrixXd inv_sigma(sigma0_.llt().solve(I));

  typedef ConditionalMode<double> TConditionalMode;
  typedef typename TConditionalMode::TVector TVector;
  typedef typename TConditionalMode::MatrixType MatrixType;

  TConditionalMode f(mu, inv_sigma);
  cppoptlib::BfgsSolver<TConditionalMode> solver;
  Rcpp::Rcout << std::endl;
  for(int i =0; i < n; i++){
    f.update_x(X_.col(i));
    Eigen::VectorXd x0_ = A_.col(i);
    solver.minimize(f, x0_);
    A_.col(i) = x0_;
  }

  Eigen::VectorXd mu_prev(A_.rowwise().mean());
  Eigen::MatrixXd centered = A_.transpose().rowwise() - A_.transpose().colwise().mean();
  Eigen::MatrixXd sigma = (centered.adjoint() * centered) / double(A_.cols() - 1);



  // Eigen::VectorXd x0(k);
  // for(int i=0; i < k; i++){
  //   x0(i) = log(x(i)+0.5) - log(x(k)+0.5);
  // }
//
//   TConditionalMode f(x_, mu_, sigma_);
//   solver.minimize(f, x0);


  // int step = 0;
  // while((max(abs(mu_prev - mu)) > tol & step < em_max_steps) | step == 0){
  //    step++;
  //    for(int i = 0; i < n; i++){
  //      A.col(i) = mvf_maximum_alr(X.col(i), mu, inv_sigma, A.col(i));
  //    }
  //    mu_prev = mu;
  //    mu = mean(A,1);
  //    inv_sigma = inv_sympd(cov(A.t()));
  // }
  return Rcpp::List::create(mu0);
}




// Rcpp::NumericVector example_NewtonDescentSolver(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericMatrix sigma) {
//   int k = mu.size();
//   Eigen::VectorXd x_(x.size());
//   Eigen::VectorXd mu_(mu.size());
//   Eigen::MatrixXd sigma_(sigma.nrow(), sigma.ncol());
//   for(int i=0; i<x.size(); i++) x_(i) = x(i);
//   for(int i=0; i<mu.size(); i++) mu_(i) = mu(i);
//   for(int i=0; i<sigma.nrow(); i++) for(int j=0; j<sigma.ncol(); j++) sigma_(i, j) = sigma(i, j);
//
//   typedef ConditionalMode<double> TConditionalMode;
//   typedef typename TConditionalMode::TVector TVector;
//   typedef typename TConditionalMode::MatrixType MatrixType;
//
//   TConditionalMode f(x_, mu_, sigma_);
//   cppoptlib::NewtonDescentSolver<TConditionalMode> solver;
//   Eigen::VectorXd x0(mu.size());
//   for(int i=0; i<mu.size(); i++){
//     x0(i) = log(x(i)+0.5) - log(x(k)+0.5);
//   }
//   solver.minimize(f, x0);
//
//   NumericVector res(x0.size());
//   for(int i=0; i<res.size(); i++) res(i) = x0(i);
//   //Rcpp::Rcout << "argmin      " << x0.transpose() << std::endl;
//   //Rcpp::Rcout << "f in argmin " << f(x0) << std::endl;
//   return res;
// }
