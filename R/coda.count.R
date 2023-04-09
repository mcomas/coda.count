#' @description
#' coda.count contains functions to estimate and work with multivariate
#' count data models. The package contains functions to deal with the
#' Dirichlet-multinomial distributions and, specially, the logratio-normal-multinomial
#' distribution (lrnm), or its particular version, the logistic-normal-multinomial.
#'
#' Different approaches to calculate or approximate the moments of a
#' lrnm are available: Montecarlo and quasi-Montecarlo
#' simulation, Gauss-Hermite quadrature and variational methods.
#'
#' To approximate the lrnm distribution the package includes a conditional version
#' were certain counts can be fixed during the estimation.
#' @useDynLib coda.count
#' @importFrom Rcpp evalCpp
"_PACKAGE"
