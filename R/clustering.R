#' @export
count.dist = function(X, mu, sigma, B = NULL){
  if(is.null(B)){
    B = coda.base::ilr_basis(ncol(X))
  }
  lN = lrnm_posterior_approx(X, mu, sigma, B)
  Bhattacharyya = function(N1, N2){
    SIGMA = (N1$sigma + N2$sigma) / 2
    invSIGMA = solve( SIGMA )
    d = 1/8 * t(N1$mu - N2$mu) %*% invSIGMA %*% (N1$mu - N2$mu) +
      1/2 * log( det(SIGMA) / sqrt(det(N1$sigma)*det(N2$sigma)) )
    d[1]
  }
  N = nrow(X)
  D = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      D[i,j] = Bhattacharyya(lN[[i]], lN[[j]])
    }
  }
  D = as.dist(D)
  D
}
