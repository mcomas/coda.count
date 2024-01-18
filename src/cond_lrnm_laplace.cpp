// // [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_DONT_PRINT_ERRORS
//
// #include "hermite.h"
// #include "dm.h"
// #include "lrnm_utils.h"
// #include "lrnm_gaussian_approx.h"
// #include "coda_base.h"
//
// using namespace Rcpp;
//
//
// // [[Rcpp::export]]
// List c_cond_lrnm_posterior_laplace(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma,
//                                    arma::mat& C){
//
//   int n = X.n_cols;
//   int D = X.n_rows;
//   int d = D - 1;
//   int dmax = std::min(d,n-1);
//
//   arma::vec d0 = arma::vec(n);
//   arma::umat I0 = arma::umat(D, n);
//   for(int k=0; k<n; k++){
//     d0(k) = arma::accu(C.col(k) == 0);
//     I0.col(k).head(d0(k)) = arma::find(C.col(k) == 0);
//   }
//
//   arma::vec d1 = d-d0;
//   arma::vec D1 = 1+d1;
//   arma::umat I1 = arma::umat(D, n);
//   for(int k=0; k<n; k++){
//     I1.col(k).head(D1(k)) = arma::find(C.col(k) > 0);
//   }
//
//   arma::cube B = arma::zeros(d,D,n);
//   for(int k=0; k<n; k++){
//     if(d0(k) > 0){
//       if(d0(k) > 1){
//         arma::uvec i0 = I0.col(k).head(d0(k));
//         for(unsigned int i = 0; i < d0(k) - 1; i++){
//           unsigned int I1 = i + 1;
//           unsigned int I2 = i + 2;
//           double l = 1/std::sqrt((double)(I1*I2));
//           double r = - std::sqrt((double)I1/I2);
//           for(unsigned int j = 0; j < I1; j++){
//             B(i,i0(j),k) = l;
//           }
//           B(i,i0(I1),k) = r;
//         }
//       }
//       for(unsigned int j = 0; j < D; j++){
//         if(C(j,k) == 0){
//           B(d0(k) - 1, j, k) = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
//         }else{
//           B(d0(k) - 1, j, k) = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
//         }
//       }
//       if(D1(k) > 1){
//         arma::uvec i1 = I1.col(k).head(D1(k));
//         for(unsigned int i = 0; i < D1(k) - 1; i++){
//           unsigned int I1 = i + 1;
//           unsigned int I2 = i + 2;
//           double l = 1/std::sqrt((double)(I1*I2));
//           double r = - std::sqrt((double)I1/I2);
//           for(unsigned int j = 0; j < I1; j++){
//             B(d0(k)+i,i1(j),k) = l;
//           }
//           B(d0(k)+i,i1(I1),k) = r;
//         }
//       }
//     }else{
//       B.slice(k) = ilr_basis_default(D1(k)).t();
//     }
//   }
//   // EASY OPTIMISATION: Construction of H1 can be incorporated when creating cube B
//   arma::mat H1 = arma::zeros(d,n);
//   for(int k=0;k<n;k++){
//     if(d1(k) > 0){
//       arma::vec Xc = X.col(k);
//       H1.col(k).head(d1(k)) = ilr_basis_default(D1(k)).t() * log(Xc(I1.col(k).head(D1(k))));
//     }
//   }
//
//   arma::vec d1s = arma::zeros(n);
//   arma::vec d0s = arma::zeros(n);
//
//   arma::cube B1s = arma::zeros(d,d,n);
//   arma::cube B0s = arma::zeros(d,d,n);
//
//   arma::mat B0s_mu = arma::zeros(d,n);
//   arma::cube B0s_sigma = arma::zeros(d,d,n);
//
//   arma::mat B0s_N_mu = arma::zeros(d,n);
//   arma::cube B0s_N_sigma = arma::zeros(d,d,n);
//
//   arma::mat clr_E1(D,n);
//   arma::cube clr_E2(D,D,n);
//
//   for(int k=0;k<n;k++){
//     arma::vec B_mu = B.slice(k) * clr_mu;
//     arma::mat B_sigma = B.slice(k) * clr_sigma * B.slice(k).t();
//     if(d0(k) > 0){
//       arma::mat B0;
//       arma::mat B0c_mu;
//       arma::mat B0c_sigma;
//       if(d1(k)==0){
//         B0c_mu = B_mu.head(d0(k));
//         B0c_sigma = B_sigma.submat(0,0,d0(k)-1,d0(k)-1);
//
//         arma::vec B0c_eigval;
//         arma::mat B0c_eigvec;
//         eig_sym(B0c_eigval, B0c_eigvec, B0c_sigma);
//         B0 = B0c_eigvec.cols(arma::find(B0c_eigval > 1e-5)).t();
//         d0s(k) = B0.n_rows;
//       }else{ // d1(k)>0
//         arma::vec h1 = H1.col(k).head(d1(k));
//
//         arma::vec B1_eigval;
//         arma::mat B1_eigvec;
//         eig_sym(B1_eigval, B1_eigvec, B_sigma.submat(d0(k), d0(k), d-1, d-1));
//         arma::mat B1 = B1_eigvec.cols(arma::find(B1_eigval > 1e-10)).t();
//         d1s(k) = B1.n_rows;
//         if(B1.n_rows >= dmax){
//           // We leave 1-dimension to fit the zero part
//           B1 = B1.tail_rows(dmax-1);
//           d1s(k) = dmax-1;
//         }
//         B1s.slice(k).submat(0,0,d1s(k)-1,d1(k)-1) = B1;
//
//         arma::mat B1_sigma = B1 * B_sigma.submat(d0(k), d0(k), d-1, d-1) * B1.t();
//         arma::mat B1_0 = B_sigma.submat(0, d0(k), d0(k)-1, d-1) * B1.t() * arma::inv_sympd(B1_sigma);
//
//         B0c_mu = B_mu.head(d0(k)) + B1_0 * B1 * (h1 - B_mu.subvec(d0(k), d-1));
//         B0c_sigma = B_sigma.submat(0,0,d0(k)-1,d0(k)-1) - B1_0 * B1 * B_sigma.submat(d0(k), 0,  d-1, d0(k)-1);
//
//         arma::vec B0c_eigval;
//         arma::mat B0c_eigvec;
//         eig_sym(B0c_eigval, B0c_eigvec, B0c_sigma);
//         B0 = B0c_eigvec.cols(arma::find(B0c_eigval > 1e-5)).t();
//         d0s(k) = B0.n_rows;
//       }
//
//       if(d0s(k)>0){
//         B0s.slice(k).submat(0, 0, d0s(k)-1, d0(k)-1) = B0;
//         arma::vec B0_mu = B0 * B0c_mu;
//         arma::mat B0_sigma = B0 * B0c_sigma * B0.t();
//
//         B0s_mu.col(k).subvec(0, d0s(k)-1) = B0_mu;
//         B0s_sigma.slice(k).submat(0, 0, d0s(k)-1, d0s(k)-1) = B0_sigma;
//
//         arma::mat inv_B0_sigma = arma::inv_sympd(B0_sigma);
//         // Finding a Gaussian approximation
//         int current_iter = 0;
//         int dh0 = d0s(k);
//         arma::vec h0(dh0);
//         arma::vec deriv1(dh0);
//         arma::mat deriv2(dh0, dh0);
//         arma::vec step = arma::ones(dh0);
//         double eps = 1e-5;
//         int max_iter = 100;
//
//         arma::mat Binv = B.slice(k).t();
//         arma::mat T0 = B0 * B.slice(k).head_rows(d0(k));
//         do{
//
//           current_iter++;
//
//           arma::vec h;
//           if(d1(k)==0){
//             h = B0.t() * h0;
//           }else{
//             arma::vec h1 = H1.col(k).head(d1(k));
//             h = arma::join_cols(B0.t() * h0, h1);
//           }
//           arma::vec eBh = exp(Binv * h);
//
//           double  w = arma::accu(eBh);
//
//           arma::vec  wBi = T0 * eBh;
//           deriv1 = -inv_B0_sigma  * (h0-B0_mu) + T0 * X.col(k) - sum(X.col(k)) * wBi / w;
//
//           arma::mat wBij(dh0,dh0);
//           for(int i=0; i<dh0; i++){
//             for(int j=0; j<dh0; j++){
//               wBij(i,j) = arma::accu(T0.row(i).t() % T0.row(j).t() % eBh);
//             }
//           }
//           deriv2 = -inv_B0_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);
//
//           step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
//           h0 = h0 - 0.9 * step;
//
//         }while( norm(step, 2) > eps && current_iter < max_iter);
//         arma::vec B0_N_mu = h0;
//         arma::mat B0_N_inv_sigma = -deriv2;
//         arma::mat B0_N_sigma = arma::inv_sympd(B0_N_inv_sigma);
//         B0s_N_mu.col(k).subvec(0, d0s(k)-1) = B0_N_mu;
//         B0s_N_sigma.slice(k).submat(0, 0, d0s(k)-1, d0s(k)-1) = B0_N_sigma;
//
//         double l_cmult = l_multinomial_const(X.col(k));
//
//         arma::vec eigval;
//         arma::mat eigvec;
//         eig_sym(eigval, eigvec, B0_N_sigma);
//         arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
//
//         unsigned int index[(unsigned)d0s(k)+1];
//         for(unsigned int i = 0; i <= d0s(k); i++) index[i] = 0;
//         int position = 0, kappa = 0;
//         arma::vec l_dens(std::pow(order, d0s(k)));
//         double l_dens_max = 0;
//
//         do{
//           double w = 0;
//           for(unsigned int i = 0; i < d0s(k); i++){
//             h0(i) = uni_hermite(index[i],0);
//             w += uni_hermite(index[i],1);
//           }
//           h0 = B0_N_mu + rotation * h0;
//           arma::vec h;
//           if(d1(k)==0){
//             h = B0.t() * h0;
//           }else{
//             arma::vec h1 = H1.col(k).head(d1(k));
//             h = arma::join_cols(B0.t() * h0, h1);
//           }
//           arma::vec Bh = B.slice(k).t() * h;
//           Bh = Bh - Bh.max();
//
//           arma::vec p = exp(Bh);
//           l_dens(kappa) = w + l_dnormal_vec(h0, B0_mu, inv_B0_sigma) -
//             l_dnormal_vec(h0, B0_N_mu, B0_N_inv_sigma) +
//             arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;
//
//           if(l_dens(kappa) > l_dens_max) l_dens_max = l_dens(kappa);
//
//
//           // Calculate next coordinate
//           index[position]++;
//           while(index[position] == order){
//             index[position] = 0;
//             position++;
//             index[position]++;
//           }
//           position = 0;
//           kappa++;
//         } while (index[(unsigned)d0s(k)] == 0);
//         arma::vec dens = exp(l_dens-l_dens_max);
//
//         double M0 = 0;
//         arma::vec M1 = arma::zeros(d0s(k));
//         arma::mat M2 = arma::zeros(d0s(k),d0s(k));
//
//         for(unsigned int i = 0; i <= d0s(k); i++) index[i] = 0;
//         position = 0;
//         kappa = 0;
//         do{
//           for(unsigned int i = 0; i < d0s(k); i++){
//             h0(i) = uni_hermite(index[i],0);
//           }
//           h0 = B0_N_mu + rotation * h0;
//           arma::vec h;
//           if(d1(k)==0){
//             h = B0.t() * h0;
//           }else{
//             arma::vec h1 = H1.col(k).head(d1(k));
//             h = arma::join_cols(B0.t() * h0, h1);
//           }
//
//           M0 += dens(kappa);
//           M1 += h0 * dens(kappa);
//           M2 += h0 * h0.t() * dens(kappa);
//
//           // Calculate next coordinate
//           index[position]++;
//           while(index[position] == order){
//             index[position] = 0;
//             position++;
//             index[position]++;
//           }
//           position = 0;
//           kappa++;
//         } while (index[(unsigned)d0s(k)] == 0);
//
//         // Expectations in the clr-space
//         Eh.subvec(0,d0(k)-1) = B0.t() * (M1/M0);
//         Ehh.submat(0,0,d0(k)-1,d0(k)-1) = B0.t() *  (M2/M0) * B0;
//         if(d1(k)>0){
//
//           arma::vec h1 = H1.col(k).head(d1(k));
//           Eh.subvec(d0(k),d-1) = h1;
//
//           Ehh.submat(0,d0(k),d0(k)-1, d-1) = Eh.subvec(0,d0(k)-1) * h1.t();
//           Ehh.submat(d0(k),0,d-1,d0(k)-1) = h1 * Eh.subvec(0,d0(k)-1).t();
//           Ehh.submat(d0(k), d0(k), d-1, d-1) = h1 * h1.t();
//
//         }
//         clr_E1.col(k) = B.slice(k).t() * Eh;
//         clr_E2.slice(k) = B.slice(k).t() * Ehh * B.slice(k);
//         //       //// HERMITE
//       }
//     }else{ // d0(k) = 0
//       arma::vec B1_eigval;
//       arma::mat B1_eigvec;
//       eig_sym(B1_eigval, B1_eigvec, clr_sigma);
//       arma::mat B1 = B1_eigvec.cols(arma::find(B1_eigval > 1e-10)).t();
//       d1s(k) = B1.n_rows;
//       clr_E1.col(k) = ilr_basis_default(D) * H1.col(k).head(d1(k));
//       clr_E2.slice(k) = clr_E1.col(k) * clr_E1.col(k).t();
//     }
//   }
//
//   // M["B0s_mu"] = B0s_mu;
//   // M["B0s_sigma"] = B0s_sigma;
//   // M["B0s_N_mu"] = B0s_N_mu;
//   // M["B0s_N_sigma"] = B0s_N_sigma;
//   // M["B0s_pc1"] = B0s_pc1;
//   // M["B0s_pc1_mu"] = B0s_pc1_mu;
//   // M["B0s_pc1_sigma"] = B0s_pc1_sigma;
//   // M["B0s_pc1_N_mu"] = B0s_pc1_N_mu;
//   // M["B0s_pc1_N_sigma"] = B0s_pc1_N_sigma;
//
//   List M;
//   M["C"] = C;
//   M["I0"] = I0;
//   M["d0"] = d0;
//   M["I1"] = I1;
//   M["d1"] = d1;
//   M["B"] = B;
//   M["H1"] = H1;
//   M["d1s"] = d1s;
//   M["d0s"] = d0s;
//   M["B1s"] = B1s;
//   M["B0s"] = B0s;
//   M["clr_E1"] = clr_E1;
//   M["clr_E2"] = clr_E2;
//   return(M);
// }
