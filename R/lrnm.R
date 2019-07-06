#' Probability mass function of a logratio-normal-multinomial distribution
#'
#' @param x count vector
#' @param mu Parameter mu of a logratio-normal-multinomial distribution
#' @param sigma Parameter sigma of a logratio-normal-multinomial distribution
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param method Methods available: 'hermite' for Hermite integration, 'mc' for Monte Carlo integration
#' @param hermite.order order of Hermite polynomials
#' @return probability mass function evaluated in x
#' @examples
#' X = apply(simplex_lattice(19,2), 2, rev)
#' sum(p <- apply(X, 1, dlrnm, -0.5, 0.1))
#' names(p) = sprintf("(%d,%d)", X[,1], X[,2])
#' barplot(p, cex.axis = 0.8, cex.names = 0.8, las=2,
#'         main = paste('Log-ratio-normal-multinomial',
#'                      'probability mass function', sep="\n"),
#'         xlab = 'Count', ylab = 'Probability')
#' @export
dlrnm = function(x, mu, sigma, B = NULL,
                 method = ifelse(length(mu) %in% 1:5, 'hermite', 'mc'),
                 hermite.order = 10){
  if(is.null(B)){
    B = coda.base::ilr_basis(length(mu)+1)
  }
  sigma = as.matrix(sigma)
  if(method == 'hermite'){
    if(is.vector(x)){
      return(c_d_lrnm_hermite(x, mu, sigma, pinv(t(B)), hermite.order))
    }else{
      return(apply(x, 1, c_d_lrnm_hermite, mu, sigma, pinv(t(B)), hermite.order))
    }
  }
  if(method == 'mc'){
    message('Method Monte Carlo not available yet')
  }
}

#' Estimate the parameters of a logratio-normal-multinomial distribution
#'
#' @param X count sample
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param probs boolean indicating if expected posterior probabilities are returned.
#' @param hermite.order order of Hermite polynomials
#' @param eps precision used for the final estimates
#' @param max_iter maximum number of iterations for the iterative procedure used to estimate the parameter
#' @return Estimated parameters mu and sigma
#' @export
fit_lrnm = function(X, B = NULL, probs = FALSE, hermite.order = 5, eps = 1e-8, max_iter = 500){
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  fit = c_fit_lm_lrnm_hermite_centered(Y = as.matrix(X), B = B, X = matrix(1, nrow(X)),
                                       order = hermite.order, eps = eps, max_iter = max_iter)
  if(fit[[4]] == max_iter){
    warning("Maximum number of iterations exhausted.")
  }
  if(probs){
    return(list('mu' = fit[[1]], 'sigma' = fit[[2]], 'P' = coda.base::composition(fit[[3]], B), iter = fit[[4]]))
  }else{
    return(list('mu' = fit[[1]], 'sigma' = fit[[2]], iter = fit[[4]]))
  }

}

#' Estimate an approximation for the posterior distribution of the log-ratio-normal-multinomial
#' distribution.
#'
#' @param X matrix of counts
#' @param mu mean parameter
#' @param sigma covariance parameter
#' @param B compositional basis
#' @param method Default 'laplace'
#' @return A cube with the parameters
#' @export
lrnm_posterior_approx = function(X, mu, sigma, B, method = 'laplace'){
  D = length(mu) + 1
  if(method == 'laplace'){
    res = c_posterior_approximation(X, mu, sigma, B)
    res = apply(res, 3, function(x) list(mu = x[,D], sigma = x[,-D]))
  }
  res
}
