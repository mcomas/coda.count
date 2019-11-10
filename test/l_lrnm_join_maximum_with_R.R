x = c(0,12,38)
mu = c(0.3790, 2.3362)
inv_sigma = matrix(c(13.0432,-3.8192,
                     -3.8192, 1.2225), 2)
Binv = ilr_basis(3)
l_lrnm_join_maximum(x, mu, inv_sigma,  Binv)

k = ncol(Binv)

if(min(x) > 0 & max(x) > 5){
  h = MASS::ginv(Binv) %*% log(x);
}else if(sum(x) > 100){
  h = MASS::ginv(Binv) %*% log(x+1);
}else{
  h = mu
}

step = rep(0, k);

current_iter = 0;
do{

  current_iter = current_iter + 1;

  deriv1 = l_lrnm_join_d1(h, x, mu, inv_sigma, Binv);
  deriv2 = l_lrnm_join_d2(h, x, mu, inv_sigma, Binv);

  Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
  step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
  h = h - 0.9 * step;
}while( norm(step, 2) > eps && current_iter < max_iter);
