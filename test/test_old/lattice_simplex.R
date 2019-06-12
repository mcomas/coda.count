Rcpp::sourceCpp('src/utils.cpp')
print("Running")
nrow(c_simplex_lattice(30,5))

