#' Probability mass function of a logratio-normal-multinomial distribution
#'
#' @param x count vector
#' @param mu Parameter mu of a logratio-normal-multinomial distribution
#' @param sigma Parameter sigma of a logratio-normal-multinomial distribution
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param method Methods available: 'hermite' for Hermite integration, 'mc' for Monte Carlo integration
#' @param hermite.order order of Hermite polynomials
#' @param mc.nsim number of low-discrepancy vectors to be generated for
#'  approximating during the E-step phase
#' @param mc.znorm matrix used to evaluate the montecarlo method, if not defined,
#'  gaussian random variables generated with R are used.
#' @return probability mass function evaluated in x
#' @examples
#' X = apply(simplex_lattice(19,2), 2, rev)
#' sum(p <- dlrnm(X, -0.5, 0.1))
#' names(p) = sprintf("(%d,%d)", X[,1], X[,2])
#' barplot(p, cex.axis = 0.8, cex.names = 0.8, las=2,
#'         main = paste('Log-ratio-normal-multinomial',
#'                      'probability mass function', sep="\n"),
#'         xlab = 'Count', ylab = 'Probability')
#' @export
dlrnm = function(x, mu, sigma, B = NULL,
                 method = ifelse(length(mu) %in% 1:4, 'hermite', 'mc'),
                 hermite.order = 10,
                 mc.nsim = 1000, mc.znorm = NULL){
  if(is.null(B)){
    B = coda.base::ilr_basis(length(mu)+1)
  }
  sigma = as.matrix(sigma)
  if(method == 'hermite'){
    if(is.vector(x)){
      return(as.vector(c_d_lrnm_hermite_mat(cbind(x), mu, sigma, pinv(t(B)), hermite.order)))
    }else{
      return(as.vector(c_d_lrnm_hermite_mat(t(x), mu, sigma, pinv(t(B)), hermite.order)))
    }
  }
  if(method == 'mc'){
    if(is.null(mc.znorm)){
      Z = t(randtoolbox::sobol(mc.nsim, dim = length(mu), normal = TRUE))
    }else{
      Z = mc.znorm
    }
    if(is.vector(x)){
      return(as.vector(c_d_lrnm_montecarlo(cbind(x), mu, sigma, pinv(t(B)), Z)))
    }else{
      return(as.vector(c_d_lrnm_montecarlo(t(x), mu, sigma, pinv(t(B)), Z)))
    }
  }
}

#' Estimate the parameters of a logratio-normal-multinomial distribution
#'
#' @param X count sample
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param method Method to use to estimate the parameters: 'hermite', 'laplace', 'mc' (default) or 'vem'
#' @param hermite.order order of Hermite polynomials
#' @param mc.nsim number of low-discrepancy vectors to be generated for
#'  approximating during the E-step phase
#' @param mc.znorm matrix used to evaluate the montecarlo method, if not defined,
#'  gaussian random variables generated with R are used.
#' @param eps precision used for the final estimates. 1e-5 for hermite method
#' and 1e-3 for montecarlo method.
#' @param max_iter maximum number of iterations for the iterative procedure
#' used to estimate the parameter
#' @param debug if set to TRUE full list with debugging parameters is returned (this list is highly possible to change in different versions)
#' @return Estimated parameters mu and sigma
#' @export
fit_lrnm = function(X, B = NULL, method = 'mc',
                    hermite.order = 5,
                    mc.nsim = 500, mc.znorm = NULL, eps = NULL, max_iter = 500,
                    B_fixed = FALSE,debug = FALSE){
  jmax = apply(X, 2, max)
  if(min(jmax) == 0){
    stop(sprintf("All observation have zero in part %d", which.min(jmax)), call. = FALSE)
  }
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  d = ncol(X)-1
  if(!method %in% c('mc','laplace','hermite', 'vem')){
    stop("Method should be 'mc', 'laplace' or 'hermite'")
  }
  if(method == 'mc'){
    if(is.null(eps)){
      eps = 0.001
    }
    if(is.null(mc.znorm)){
      Z = randtoolbox::sobol(mc.nsim, dim = ncol(X)-1, normal = TRUE)
    }else{
      Z = mc.znorm
    }
    if(B_fixed){
      fit = c_lrnm_fit_fixed_montecarlo(t(X), t(Z), eps, max_iter)
    }else{
      fit = c_lrnm_fit_montecarlo(t(X), t(Z), eps, max_iter)
    }

  }
  if(method=='laplace'){
    if(is.null(eps)){
      eps = 1e-05
    }
    fit = c_lrnm_fit_laplace(t(X), eps, max_iter)
  }
  if(method=='hermite'){
    if(is.null(eps)){
      eps = 1e-05
    }
    if(B_fixed){
      fit = c_lrnm_fit_fixed_hermite(t(X), hermite.order, eps, max_iter)
    }else{
      fit = c_lrnm_fit_hermite(t(X), hermite.order, eps, max_iter)
    }

  }
  if(method=='vem'){
    if(is.null(eps)){
      eps = 1e-05
    }
    fit = c_vem_lrnm_fit(t(X), eps, max_iter)
  }
  if(debug) return(fit)
  result = list('mu' = (t(B) %*% fit$clr_mu)[,1], 'sigma' = t(B) %*% fit$clr_sigma %*% B,
                'P' = coda.base::composition(t(fit$clr_E1), 'clr'),
                'iter' = fit$em_iter, eps = eps)

  if(result$iter == max_iter){
    warning("Maximum number of iterations exhausted.")
  }
  return(result)

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
    res = apply(X, 1, c_lrnm_posterior_approximation_vec, mu, sigma, B, eps = 0.001, niter = 100,
          simplify = FALSE)
    res = lapply(res, function(x) list(mu = x[,D], sigma = x[,-D]))
  }
  if(method == 'vem'){
    stop("Not implemented yet")
  }
  res
}

