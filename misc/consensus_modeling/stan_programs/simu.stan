transformed data {
  int<lower=0> K = 3;
  int<lower=0> N_reviewers = 5;
  int<lower=0> N_assessments = 110;
  
  vector[K] psi = [ 0.4, 0.3, 0.3 ]';
  vector[K] fidelity[N_reviewers, K];
  
  for (r in 1:N_reviewers) {
    for (k in 1:K) {
      vector[K] alpha = rep_vector(1, K);
      alpha[k] = 10;
      fidelity[r, k] = dirichlet_rng(alpha);
    }
  }
}

generated quantities {
  int<lower=0, upper=K> z[N_assessments];
  int<lower=0, upper=K> y[N_assessments, N_reviewers];
  
  for (n in 1:N_assessments) {
    z[n] = categorical_rng(psi);
    
    for (r in 1:N_reviewers) {
      y[n, r] = categorical_rng(fidelity[r, z[n]]);
    }
  }
}
