#' Probability mass function of a logratio-normal-multinomial distribution
#'
#' @param x count vector
#' @param mu Parameter mu of a logratio-normal-multinomial distribution
#' @param sigma Parameter sigma of a logratio-normal-multinomial distribution
#' @param method Methods available: 'hermite' for Hermite integration, 'mc' for Monte Carlo integration
#' @return probability mass function evaluated in x
#' @examples
#' X = apply(simplex_lattice(10,2), 2, rev)
#' sum(p <- apply(X, 1, dlrnm, -0.5, 0.1))
#' names(p) = sprintf("(%d,%d)", X[,1], X[,2])
#' barplot(p, cex.axis = 0.8, cex.names = 0.8)
#' @export
dlrnm = function(x, mu, sigma, method = ifelse(length(mu) %in% 1:2, 'hermite', 'mc'),
                 hermite.order = 100, hermite.step_by = 100, hermite.eps = 1e-06,
                 hermite.max_steps = 10){
  sigma = as.matrix(sigma)
  if(method == 'hermite'){
    if(is.vector(x)){
      return(c_dlrnm_hermite(x, mu, sigma, hermite.order, hermite.step_by, hermite.eps, hermite.max_steps))
    }else{
      return(apply(X, 1, c_dlrnm_hermite, mu, sigma, hermite.order, hermite.step_by, hermite.eps, hermite.max_steps))
    }
  }
  if(method == 'mc'){
    message('Method Monte Carlo not available yet')
  }
}

#' Estimate the parameters of a logratio-normal-multinomial distribution
#'
#' @param X count sample
#' @return Estimated parameters mu and sigma
#' @export
fit_lrnm = function(X, eps = 0.0001, maxiter = 5000){
  message('Method not available yet')
}

