#' Estimate the parameters of a logratio-normal-multinomial distribution
#' when some logratios are fixed.
#'
#' @param X count sample
#' @param C conditioning variable TRUE variables are fixed
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param method Method to use to estimate the parameters: 'hermite', 'laplace', 'mc' (default) or 'vem'
#' @param low.dim Number of dimensions used to calculate expected values, 0: all dimensions.
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
fit_conditional_lrnm = function(X, C = X > 0, B = NULL, method = 'mc',
                                low.dim = 0, hermite.order = 5,
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
  if(!method %in% c('mc','laplace','hermite')){
    stop("Method should be 'mc', 'laplace' or 'hermite'")
  }
  if(low.dim == 0){
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
        fit = c_cond_lrnm_fit_fixed_montecarlo(t(X), t(C), t(Z), eps, max_iter)
      }else{
        fit = c_cond_lrnm_fit_montecarlo(t(X), t(C), t(Z), eps, max_iter)
      }


    }
    if(method=='laplace'){
      if(is.null(eps)){
        eps = 1e-05
      }
      message("Not implemented yet")
      # fit = c_lrnm_fit_laplace(t(X), eps, max_iter)
    }
    if(method=='hermite'){
      if(is.null(eps)){
        eps = 1e-05
      }
      if(B_fixed){
        fit = c_cond_lrnm_fit_fixed_hermite(t(X), t(C), hermite.order, eps, max_iter)
      }else{
        fit = c_cond_lrnm_fit_hermite(t(X), t(C), hermite.order, eps, max_iter)
      }

    }
  }else{
    if(is.null(eps)){
      eps = 1e-05
    }
    fit = c_low_dim_cond_lrnm_fit_hermite(t(X), t(C), low.dim, hermite.order, eps, max_iter)
  }
  if(debug) return(fit)
  result = list('mu' = (t(B) %*% fit$clr_mu)[,1], 'sigma' = t(B) %*% fit$clr_sigma %*% B,
                'P' = coda.base::composition(t(fit$clr_H), 'clr'),
                'iter' = fit$em_iter, eps = eps)


  if(result$iter == max_iter){
    warning("Maximum number of iterations exhausted.")
  }
  return(result)
}

#' @export
cond_lrnm_posterior_moments_hermite = c_cond_lrnm_posterior_moments_hermite
