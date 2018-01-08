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
