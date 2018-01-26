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

// [[Rcpp::export]]
arma::mat c_lrnm_fit_maximum_alr(arma::mat X, arma::vec mu0, arma::mat sigma0,
                                  double tol = 10e-6, int em_max_steps = 100,
                                  double min_evalue = 0.0001){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  Eigen::MatrixXd X_(K, n);
  for(int i = 0; i < K; i++) for(int j = 0; j < n; j++) X_(i,j) = X(j,i);

  Eigen::VectorXd mu0_(k);
  for(int i = 0; i < k; i++) mu0_(i) = mu0(i);

  Eigen::MatrixXd sigma0_(k, k);
  for(int i = 0; i < k; i++) for(int j = 0; j < k; j++) sigma0_(i, j) = sigma0(i, j);

  Eigen::MatrixXd A_(k, n);
  for(int i = 0; i < n; i++){
    double denom = log(X_(k,i) + 0.5);
    for(int j = 0; j < k; j++){
      A_(j,i) = log(X_(j,i) + 0.5) - denom;
    }
  }
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
  Eigen::VectorXd mu(mu0_);
  Eigen::MatrixXd inv_sigma(sigma0_.llt().solve(I));

  typedef ConditionalMode<double> TConditionalMode;
  typedef typename TConditionalMode::TVector TVector;
  typedef typename TConditionalMode::MatrixType MatrixType;

  TConditionalMode f0(mu, inv_sigma);
  cppoptlib::BfgsSolver<TConditionalMode> solver;
  for(int i =0; i < n; i++){
    f0.update_x(X_.col(i));
    Eigen::VectorXd x0_ = A_.col(i);
    solver.minimize(f0, x0_);
    A_.col(i) = x0_;
  }

  Eigen::VectorXd mu_prev = mu.replicate(1, 1);
  mu = A_.rowwise().mean();
  Eigen::MatrixXd centered = A_.transpose().rowwise() - A_.transpose().colwise().mean();
  Eigen::MatrixXd sigma_ini = (centered.adjoint() * centered) / double(A_.cols() - 1);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(sigma_ini);
  Eigen::MatrixXd D = es.eigenvalues().array().max(min_evalue).matrix().asDiagonal();
  Eigen::MatrixXd V = es.eigenvectors();
  Eigen::MatrixXd sigma = V * D * V.transpose();

  int step = 0;
  while( ((mu_prev - mu).lpNorm<Eigen::Infinity>() > tol & step < em_max_steps) | step == 0){
    //Rcpp::Rcout << A_ << std::endl << sigma << std::endl << mu << std::endl;
    inv_sigma = sigma.llt().solve(I);
    TConditionalMode f(mu, inv_sigma);
    step++;
    for(int i =0; i < n; i++){
      f.update_x(X_.col(i));
      Eigen::VectorXd x0_ = A_.col(i);
      //Rcpp::Rcout << mu << std::endl << inv_sigma << std::endl << x0_ << std::endl;
      solver.minimize(f, x0_);
      //Rcpp::Rcout << x0_ <<std::endl <<std::endl;
      A_.col(i) = x0_;
    }
    //Rcpp::Rcout << A_ << std::endl << inv_sigma << std::endl;
    mu_prev = mu.replicate(1, 1);
    mu = A_.rowwise().mean();
    Eigen::MatrixXd centered = A_.transpose().rowwise() - A_.transpose().colwise().mean();
    sigma_ini = (centered.adjoint() * centered) / double(A_.cols() - 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(sigma_ini);
    D = es.eigenvalues().array().max(min_evalue).matrix().asDiagonal();
    V = es.eigenvectors();
    sigma = V * D * V.transpose();
  }

  arma::mat A(n, k);
  for(int i=0; i<n; i++){
    for(int j=0; j<k; j++){
      A(i,j) = A_(j,i);
    }
  }
  return A;
}

// [[Rcpp::export]]
arma::mat c_lrnm_fit_maximum_alr_centered(arma::mat X, arma::vec mu0, arma::mat sigma0,
                                          double tol = 10e-6, int em_max_steps = 100,
                                          double min_evalue = 0.0001){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  Eigen::MatrixXd X_(K, n);
  for(int i = 0; i < K; i++) for(int j = 0; j < n; j++) X_(i,j) = X(j,i);

  Eigen::VectorXd mu0_(k);
  for(int i = 0; i < k; i++) mu0_(i) = mu0(i);

  Eigen::MatrixXd sigma0_(k, k);
  for(int i = 0; i < k; i++) for(int j = 0; j < k; j++) sigma0_(i, j) = sigma0(i, j);

  Eigen::MatrixXd A_(k, n);
  for(int i = 0; i < n; i++){
    double denom = log(X_(k,i) + 0.5);
    for(int j = 0; j < k; j++){
      A_(j,i) = log(X_(j,i) + 0.5) - denom;
    }
  }
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
  Eigen::VectorXd mu(mu0_);
  Eigen::MatrixXd inv_sigma(sigma0_.llt().solve(I));

  typedef ConditionalMode<double> TConditionalMode;
  typedef typename TConditionalMode::TVector TVector;
  typedef typename TConditionalMode::MatrixType MatrixType;

  TConditionalMode f0(mu, inv_sigma);
  cppoptlib::BfgsSolver<TConditionalMode> solver;
  for(int i =0; i < n; i++){
    f0.update_x(X_.col(i));
    Eigen::VectorXd x0_ = A_.col(i);
    solver.minimize(f0, x0_);
    A_.col(i) = x0_;
  }

  Eigen::VectorXd mu_prev = mu.replicate(1, 1);
  mu = A_.rowwise().mean();
  Eigen::MatrixXd centered = A_.transpose().rowwise() - A_.transpose().colwise().mean();
  Eigen::MatrixXd sigma_ini = (centered.adjoint() * centered) / double(A_.cols() - 1);
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(k, k) * std::max(min_evalue, sigma_ini.diagonal().mean());

  int step = 0;
  while( ((mu_prev - mu).lpNorm<Eigen::Infinity>() > tol & step < em_max_steps) | step == 0){
    //Rcpp::Rcout << A_ << std::endl << sigma << std::endl << mu << std::endl;
    inv_sigma = sigma.llt().solve(I);
    TConditionalMode f(mu, inv_sigma);
    step++;
    for(int i =0; i < n; i++){
      f.update_x(X_.col(i));
      Eigen::VectorXd x0_ = A_.col(i);
      //Rcpp::Rcout << mu << std::endl << inv_sigma << std::endl << x0_ << std::endl;
      solver.minimize(f, x0_);
      //Rcpp::Rcout << x0_ <<std::endl <<std::endl;
      A_.col(i) = x0_;
    }
    //Rcpp::Rcout << A_ << std::endl << inv_sigma << std::endl;
    mu_prev = mu.replicate(1, 1);
    mu = A_.rowwise().mean();
    Eigen::MatrixXd centered = A_.transpose().rowwise() - A_.transpose().colwise().mean();
    sigma_ini = (centered.adjoint() * centered) / double(A_.cols() - 1);
    sigma = Eigen::MatrixXd::Identity(k, k) * std::max(min_evalue, sigma_ini.diagonal().mean());
  }

  arma::mat A(n, k);
  for(int i=0; i<n; i++){
    for(int j=0; j<k; j++){
      A(i,j) = A_(j,i);
    }
  }
  return A;
}
