#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int lattice_elements(int K, int SIZE){
  return R::choose(K+SIZE-1,SIZE);
}

//' @export
// [[Rcpp::export]]
NumericMatrix lattice_simplex(int K, int SIZE) {
  int nrow = lattice_elements(K,SIZE);
  int k = K - 1;
  int k1 = k - 1;
  NumericMatrix out(nrow, K);
  NumericVector x(K);
  x(0) = SIZE;

  int target = 0;
  int i = -1;
  do{
    i++;
    out(i,_) = x;
    x(target) = x(target) - 1;
    if(target < k1){
      target = target + 1;
      x(target) = 1 + x(k);
      x(k) = 0;
    }else{
      x(k) = x(k) + 1;
      while(x(target) == 0){
        target = target - 1;
        if(target == -1){
          i++;
          out(i,_) = x;
          return(out);
        }
      }
    }
  }while(true);
  return(out);
}




// function (n, m, tol = 1e-08)
// {
//   len <- max(length(n), length(m))
//   out <- numeric(len)
//   n <- rep(n, length = len)
//   m <- rep(m, length = len)
//   mint <- (trunc(m) == m)
//   out[!mint] <- NA
//   out[m == 0] <- 1
//   whichm <- (mint & m > 0)
//   whichn <- (n < 0)
//   which <- (whichm & whichn)
//   if (any(which)) {
//     nnow <- n[which]
//     mnow <- m[which]
//     out[which] <- ((-1)^mnow) * Recall(mnow - nnow - 1, mnow)
//   }
//   whichn <- (n > 0)
//     nint <- (trunc(n) == n)
//     which <- (whichm & whichn & !nint & n < m)
//     if (any(which)) {
//       nnow <- n[which]
//       mnow <- m[which]
//       foo <- function(j, nn, mm) {
//         n <- nn[j]
//         m <- mm[j]
//         iseq <- seq(n - m + 1, n)
//         negs <- sum(iseq < 0)
//         ((-1)^negs) * exp(sum(log(abs(iseq))) - lgamma(m +
//           1))
//       }
//       out[which] <- unlist(lapply(seq(along = nnow), foo, nn = nnow,
//                                   mm = mnow))
//     }
//     which <- (whichm & whichn & n >= m)
//       nnow <- n[which]
//     mnow <- m[which]
//     out[which] <- exp(lgamma(nnow + 1) - lgamma(mnow + 1) - lgamma(nnow -
//       mnow + 1))
//       nna <- !is.na(out)
//       outnow <- out[nna]
//     rout <- round(outnow)
//       smalldif <- abs(rout - outnow) < tol
//       outnow[smalldif] <- rout[smalldif]
//     out[nna] <- outnow
//       out
