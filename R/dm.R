#' Probability mass function of a Dirichlet-multinomial distribution
#'
#' @param x count vector
#' @param alpha Parameter alpha of a Dirichlet-Multinomial distribution
#' @return probability mass function evaluated in x
#' @examples
#' X = apply(simplex_lattice(19,2), 2, rev)
#' sum(p <- apply(X, 1, ddm, c(1.5,2)))
#' names(p) = sprintf("(%d,%d)", X[,1], X[,2])
#' barplot(p, cex.axis = 0.8, cex.names = 0.8, las=2,
#'         main = paste('Dirichlet-multinomial',
#'                      'probability mass function', sep="\n"),
#'                      xlab = 'Count', ylab = 'Probability')
#' D = 3
#' X = simplex_lattice(10,3)
#' (pmf <- cbind(X, ddm(X, c(0.5,0.5,0.5))))
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
    apply(x, 1, c_ddm, alpha)
  }

}

#' Estimate the parameters of a Dirichlet-multinomial distribution
#'
#' @param X count sample
#' @param eps precision used for the final estimates
#' @param max_iter maximum number of iterations for the iterative procedure used to estimate the parameter
#' @return Estimated parameters alpha
#' @examples
#' set.seed(1)
#' X = rdm(n = 1000, size = 100, c(1,4,0.5))
#' fit_dm(X)
#' @export
fit_dm = function(X, eps = 0.0001, max_iter = 500){
  fitting = c_dm_fit(t(X), eps, max_iter)
  if(fitting[["iter"]] >= max_iter){
    warning("Maximum number of iterations reached. Try to increase max_iter.")
  }
  fitting[["alpha"]][,1]
}

