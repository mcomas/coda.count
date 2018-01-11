#' Probability mass function of a Dirichlet-multinomial distribution
#'
#' @param x count vector
#' @param alpha Parameter alpha of a Dirichlet-Multinomial distribution
#' @return probability mass function evaluated in x
#' @examples
#' X = lattice_simplex(3,10)
#' (pmf <- cbind(X, apply(X, 1, ddm, c(0.5,0.5,0.5))))
#' sum(pmf[,4])
#' @export
ddm <- function(x, alpha) {
  x = trunc(x)
  dimension = ifelse(is.vector(x), length(x), ncol(x))
  if( dimension != length(alpha) ){
    stop("count and parameters must have same length")
  }
  if( ! length(alpha) > 1){
    stop("Dimension must be at least two")
  }
  if( any(x < 0) ){
    stop("All components must be non-negative")
  }
  if( !all(alpha > 0) ){
    stop("All parameters must be positive")
  }
  if(is.vector(x)){
    c_ddm(x, alpha)
  }else{
    apply(X, 1, c_ddm, alpha)
  }

}

#' Estimate the parameters of a Dirichlet-multinomial distribution
#'
#' @param X count sample
#' @return Estimated parameters alpha
#' @examples
#' set.seed(1)
#' X = rdm(c(1,4,0.5), 100, n = 100)
#' fit_dm(X)
#' @export
fit_dm = function(X, eps = 0.0001, maxiter = 5000){
  fitting = c_dm_fit(X, eps, maxiter)
  alpha = fitting[[1]]
  attr(alpha, 'iter') = fitting[[2]]
  alpha
}

