#' Multinomial random sample
#'
#' @param n number of random samples. Ignored if p is a matrix
#' @param size vector to set the multinomial sampling size. vector is reused to have length equal to the number of rows of probs
#' @param p vector or matrix specifiying the probability of each class. If a matrix is provided the number of rows
#' will be equal to the number of samples to generate. Row number i defines the multinomial probability
#' to be used when generating sample i. In this last scenario, parameter names (size and p) should be provided when calling the function.
#' @return multinomial random sample
#' @examples
#' rmultinomial(10, 4, c(0.3, 0.4, 0.3))
#' probs = matrix(c(0.2, 0.4, 0.4,
#'                  0.5, 0.1, 0.4), byrow = TRUE, nrow = 6, ncol = 3)
#' size = c(10, 500, 10000)
#' rmultinomial(size = size, p = probs)
#' @export
rmultinomial = function(n = NULL, size, p){
  if(is.matrix(p)){
    P = p
  }
  if(is.vector(p)){
    ncols_ = length(p)
    nrows_ = max(n, length(size))
    P = matrix(p, ncol = ncols_, nrow = nrows_, byrow = TRUE)
  }
  N = rep(size, length.out = nrow(P))

  c_rmultinomial_Rcpp(P, N)
}

#' Dirichlet-multinomial random sample
#'
#' @param alpha Dirichlet-multinomial parameter
#' @param size vector to set the multinomial sampling size. vector is reused to have length equal parameter n
#' @param n sample size
#' @param probs logical indicating whether multinomial probabilities should be returned
#' @return Dirichlet-multinomial random sample
#' @examples
#' rdm(100, 1000, c(1,1,1))
#' rdm(size = c(1000, 100, 10, 2, 1), alpha = c(1,1,1))
#' @export
rdm = function(n = NULL, size, alpha, probs = FALSE){
  if(is.null(n)){
    n = length(size)
  }
  N = rep(size, length.out = n)
  gen = c_rdirichletmultinomial(alpha, N)
  X = gen[[2]]
  if(probs){
    attr(X, 'probs') = gen[[1]]
  }
  X
}

#' logratio-normal-multinomial random sample
#'
#' @param n sample size. if ommited, the value is obtained as the length of vector size (in that case, all parameters should be called using their name).
#' @param size vector to set the multinomial sampling size. vector is reused to have length equal parameter n
#' @param mu logratio-normal-multinomial mu parameter
#' @param sigma logratio-normal-multinomial mu parameter
#' @param B basis in which the parameters are expressed with
#' @param probs logical indicating whether multinomial probabilities should be returned
#' @return logratio-normal-multinomial random sample
#' @examples
#' mu = c(0,1)
#' sigma = diag(2)
#' rlrnm(10, 1000, mu, sigma)
#' @export
rlrnm = function(n = NULL, size, mu, sigma, B = NULL, probs = FALSE){
  if(is.null(B)){
    B = coda.base::ilr_basis(length(mu)+1)
  }
  if(is.null(n)){
    n = length(size)
  }
  N = rep(size, length.out = n)
  gen = c_rnormalmultinomial(mu, sigma, N, pinv(B))
  X = gen[[2]]
  if(probs){
    attr(X, 'probs') = gen[[1]]
  }
  X
}

#' logratio-normal random sample defined on the Simplex
#'
#' @param n sample size
#' @param mu logratio-normal-multinomial mu parameter
#' @param sigma logratio-normal-multinomial mu parameter
#' @param B basis in which the parameters are expressed with
#' @return logratio-normal-multinomial random sample
#' @examples
#' mu = c(0,1)
#' sigma = diag(2)
#' P.ilr <- rlrnormal(100, mu, sigma)
#' colMeans(coda.base::coordinates(P.ilr))
#' cov(coda.base::coordinates(P.ilr))
#'
#' B = alr_basis(3)
#' P.alr <- rlrnormal(100, mu, sigma, B)
#' colMeans(coda.base::coordinates(P.alr, B))
#' cov(coda.base::coordinates(P.alr, B))
#' @export
rlrnormal = function(n, mu, sigma, B = NULL){
  if(is.null(B)){
    B = coda.base::ilr_basis(length(mu)+1)
  }
  c_rnormalSimplex(n, mu, sigma, pinv(B))
}
