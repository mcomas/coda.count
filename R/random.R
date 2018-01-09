#' Multinomial random sample
#'
#' @param probs matrix with number of rows equal to the number of samples to generate. Row number i contains the multinomial probability
#' to be  used when generating sample i.
#' @param size vector to set the multinomial sampling size. vector is reused to have length equal to the number of rows of probs
#' @return multinomial random sample
#' @examples
#' probs = matrix(c(0.2, 0.4, 0.4,
#'                  0.5, 0.1, 0.4), byrow = TRUE, nrow = 6, ncol = 3)
#' size = c(10, 500, 10000)
#' rmultinomial(probs, size)
#' @export
rmultinomial = function(probs, size){
  dimprobs = dim(probs)
  if( length(dimprobs) != 2 ){
    stop("probs must have dimension two (rows and columns)")
  }
  P = as.matrix(probs)
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
#' rdm(c(1,1,1), 1000, n = 100)
#' @export
rdm = function(alpha, size, n = length(size), probs = FALSE){
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
#' @param mu logratio-normal-multinomial mu parameter
#' @param sigma logratio-normal-multinomial mu parameter
#' @param size vector to set the multinomial sampling size. vector is reused to have length equal parameter n
#' @param n sample size
#' @param probs logical indicating whether multinomial probabilities should be returned
#' @return logratio-normal-multinomial random sample
#' @examples
#' mu = c(0,1)
#' sigma = diag(2)
#' rlrnm(mu, sigma, 1000, 10)
#' @export
rlrnm = function(mu, sigma, size, n = length(size), probs = FALSE, B = NULL){
  N = rep(size, length.out = n)
  if(is.null(B)){
    mu.ilr = mu
    sigma.ilr = sigma
  }else{
    B0 = ilr_basis(length(mu)+1)
    mu.ilr = mu %*% MASS::ginv(B) %*% B0
    sigma.ilr = t(B0) %*% t(MASS::ginv(B)) %*% sigma  %*% MASS::ginv(B) %*% B0
  }
  gen = c_rnormalmultinomial(mu.ilr, sigma.ilr, N)
  X = gen[[2]]
  if(probs){
    attr(X, 'probs') = gen[[1]]
  }
  X
}
