#' Create all elements of an (n,D)-simplex lattice
#'
#' @param n total sum
#' @param D simplex dimension
#' @return matrix containing all elements of an (n,D)-simplex lattice
#' @examples
#' simplex_lattice(4,3)
#' @export
simplex_lattice = function(n, D){
  c_simplex_lattice(D, n)
}
