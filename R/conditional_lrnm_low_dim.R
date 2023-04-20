#' Estimate the parameters of a logratio-normal-multinomial distribution
#' when some logratios are fixed.
#'
#' @param X count sample
#' @param C conditioning variable TRUE variables are fixed
#' @param V a matrix (same size as X) containing the directions (log-contrast) used to find the expected values in the E-step for each observation.
#' @param V.direction  'center' if sub-variety passing through the normal's center (default) or 'equal' using the equal weight for part.
#' @param B clr-basis in which mu and sigma are interpreted with. Default basis is given by coda.base::ilr_basis(length(x))
#' @param hermite.order order of Hermite polynomials
#' @param eps precision used for the final estimates. 1e-5 for Hermite method
#' and 1e-3 for Montecarlo method.
#' @param max_iter maximum number of iterations for the iterative procedure
#' used to estimate the parameter
#' @param debug if set to TRUE full list with debugging parameters is returned (this list is highly possible to change in different versions)
#' @return Estimated parameters mu and sigma
#' @export
fit_one_dimensional_conditional_lrnm = function(X, C = X > 0, V = NULL, V.direction = "center", B = NULL,
                                                hermite.order = 5, eps = 1e-5, max_iter = 500,
                                                debug = FALSE){
  jmax = apply(X, 2, max)
  if(min(jmax) == 0){
    stop(sprintf("All observation have zero in part %d", which.min(jmax)), call. = FALSE)
  }
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  d = ncol(X)-1
  if(is.null(V) & !V.direction %in% c('center','equal')){
    stop("Direction V should be specified or V.direction set to 'center' or 'equal'")
  }
  if(is.null(V)){
    if(V.direction == 'equal'){
      sel = rowSums(X==0)>0
      PART = 2*t(X==0)-1
      Vclr = matrix(0, nrow = ncol(X), ncol = nrow(X))
      Vclr[,sel] = coda.base::sbp_basis(PART[,sel], silent = TRUE)
      fit = c_cond_lrnm_V_fit_hermite(t(X), t(C), Vclr, hermite.order, eps, max_iter)
    }else{
      fit = c_cond_lrnm_one_dim_fit_hermite(t(X), t(C), hermite.order, eps, max_iter)
    }
  }else{
    fit = c_cond_lrnm_V_fit_hermite(t(X), t(C), t(V.direction), hermite.order, eps, max_iter)
  }
  if(debug) return(fit)
  result = list('mu' = (t(B) %*% fit$clr_mu)[,1], 'sigma' = t(B) %*% fit$clr_sigma %*% B,
                'P' = coda.base::composition(t(fit$clr_E1), 'clr'),
                'iter' = fit$em_iter, eps = eps)

}
