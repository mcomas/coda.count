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
#' @param probs boolean indicating if expected posterior probabilities are returned.
#' @param method Method to use to estimate the parameters: 'hermite', 'laplace', 'montecarlo' (default)
#' @param H.ini Coordinates used to initialise the expected posterior compositions given count data
#' @param mc.nsim number of samples in the Montecarlo integration process.
#' @param hermite.order order of Hermite polynomials
#' @param mc.znorm matrix used to evaluate the montecarlo method, if not defined, gaussian random variables generated with R are used.
#' @param eps precision used for the final estimates. 1e-5 for hermite method and 1e-3 for montecarlo method.
#' @param max_iter maximum number of iterations for the iterative procedure used to estimate the parameter
#' @return Estimated parameters mu and sigma
#' @export
fit_lrnm = function(X, B = NULL, probs = FALSE, method = 'mc',
                    hermite.order = 5,
                    mc.nsim = 500, mc.znorm = NULL, eps = NULL, max_iter = 500){
  jmax = apply(X, 2, max)
  if(min(jmax) == 0){
    stop(sprintf("All observation have zero in part %d", which.min(jmax)), call. = FALSE)
  }
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  d = ncol(X)-1
  if(!method %in% c('mc','laplace','hermite')){
    stop("Method should be 'mc', 'laplace' or 'hermite'")
  }
  if(method == 'mc'){
    if(is.null(eps)){
      eps = 0.001
    }
    if(is.null(mc.znorm)){
      Z = randtoolbox::sobol(mc.nsim, dim = length(mu), normal = TRUE)
    }else{
      Z = mc.znorm
    }
    fit = c_lrnm_fit_montecarlo(t(X), t(Z), eps, max_iter)

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
    fit = c_lrnm_fit_hermite(t(X), hermite.order, eps, max_iter)
  }
  if(probs){
    result = list('mu' = (t(B) %*% fit$clr_mu)[,1], 'sigma' = t(B) %*% fit$clr_sigma %*% B,
                  'P' = coda.base::composition(t(fit$clr_E1), 'clr'),
                  'iter' = fit$em_iter, eps = eps)
  }else{
    result = list('mu' = (t(B) %*% fit$clr_mu)[,1], 'sigma' = t(B) %*% fit$clr_sigma %*% B,
                  'P' = coda.base::composition(t(fit$clr_E1), 'clr'),
                  'iter' = fit$em_iter, eps = eps)
  }

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
  if(method == 'variational'){
    res = list()
    for(i in 1:nrow(X)){
      x = X[i,]
      m = coda.base::coordinates(0.5 * coda.base::composition(mu, B) + 0.5*x/sum(x), B)
      V = rep(1, length(m))
      xi = 1
      EPS = 1
      ITER = 1
      while(ITER < 50 & EPS > 0.0001){
        opt_pars = optimise_xi_m_V(m, V, xi, x, mu, sigma, B)

        ITER = ITER + 1
        EPS = max(m - opt_pars$m)^2

        xi = opt_pars$xi
        m = opt_pars$m
        V = opt_pars$V

      }
      res[[i]] = list(mu = m, sigma = V)
    }
  }
  res
}

