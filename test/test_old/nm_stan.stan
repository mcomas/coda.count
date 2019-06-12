functions {
  real mvn_lpdf(matrix h, matrix x, matrix mu, matrix sigma, matrix p, real lm_const){
    return lm_const +
           sum(x .* log(p)) +
           (- 0.5 * log_determinant(2 * pi() * sigma)
            - 0.5 * (h-mu)' * inverse(sigma) * (h-mu))[1,1];
  }
}
data {
  int K;
  matrix[K,1] x;
  matrix[K-1,1] mu;
  matrix[K-1,K-1] sigma;
  matrix[K,K-1] B;
  real MCONSTANT;
}
parameters {
  matrix[K-1,1] h;
}
transformed parameters {
  matrix[K,1] p;
  p = exp(B * h);
  p = p / sum(p);
}
model {
  h ~ mvn_lpdf(x, mu, sigma, p, MCONSTANT);
}
