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

#' @export
log_join_lrnm = function(x, h, mu, sigma, B = NULL, constant = TRUE){
  if(constant){
    l_lrnm_join_vec(h, x, mu, solve(sigma), pinv(t(B)))
  }else{
    l_lrnm_join_no_constant_vec(h, x, mu, solve(sigma), pinv(t(B)))
  }
}

#' Estimate the parameters of a logratio-normal-multinomial distribution
#'
#' @param X count sample
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param probs boolean indicating if expected posterior probabilities are returned.
#' @param montecarlo.n number of samples in the Montecarlo integration process.
#' @param hermite.order order of Hermite polynomials
#' @param eps precision used for the final estimates. 1e-5 for hermite method and 1e-3 for montecarlo method.
#' @param max_iter maximum number of iterations for the iterative procedure used to estimate the parameter
#' @return Estimated parameters mu and sigma
#' @export
fit_lrnm = function(X, B = NULL, probs = FALSE, method = 'montecarlo', H.ini = NULL,
                    montecarlo.n = 500, hermite.order = 5, Z = NULL, eps = NULL, max_iter = 500){
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  if(is.null(H.ini)){
    alpha = fit_dm(X)[,1]
    H.ini = log(t(t(as.matrix(X)) + alpha)) %*% B
  }
  d = ncol(X)-1
  if(method=='hermite'){
    if(is.null(eps)){
      eps = 1e-05
    }
    fit = c_fit_lrnm_lm_hermite(Y = as.matrix(X), B = B, X = matrix(1, nrow(X)),
                                order = hermite.order, eps = eps, max_iter = max_iter, H0 = H.ini)
  }
  if(method=='laplace'){
    if(is.null(eps)){
      eps = 1e-05
    }
    fit = c_fit_lrnm_lm_laplace(Y = as.matrix(X), B = B, X = matrix(1, nrow(X)),
                                eps = eps, max_iter = max_iter, H0 = H.ini)
  }
  if(method=='montecarlo'){
    if(is.null(eps)){
      eps = 1e-3
    }
    if(is.null(Z)){
      Z = matrix(rnorm(montecarlo.n*d), ncol = d)
      Z = rbind(Z,-Z)
    }
    fit = c_fit_lrnm_lm_montecarlo(Y = as.matrix(X), B = B, X = matrix(1, nrow(X)),
                                   Z = Z, eps = eps, max_iter = max_iter, H0 = H.ini)
  }
  if(fit[[4]] == max_iter){
    warning("Maximum number of iterations exhausted.")
  }
  if(probs){
    return(list('mu' = fit[[1]], 'sigma' = fit[[2]], 'P' = coda.base::composition(fit[[3]], B),
                iter = fit[[4]], eps = eps))
  }else{
    return(list('mu' = fit[[1]], 'sigma' = fit[[2]],
                iter = fit[[4]], eps = eps))
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
  if(method == 'variational'){
    res = list()
    for(i in 1:nrow(X)){
      x = X[i,]
      m = coordinates(0.5 * composition(mu, B) + 0.5*x/sum(x), B)
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

