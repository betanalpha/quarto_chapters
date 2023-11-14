data {
  int<lower=0> K;
  int<lower=0> N_reviewers;
  int<lower=0> N_assessments;

  int<lower=0, upper=K> y[N_assessments, N_reviewers];
  
  int<lower=0> N_calib_assessments;
  int<lower=0, upper=K> y_calib[N_calib_assessments, N_reviewers];
  int<lower=0, upper=K> z_calib[N_calib_assessments];
}

parameters {
  // One simplex for the true answers of every assessment
  simplex[K] psi[N_assessments];
  
  // One simplex per reviewer per possible true answers
  simplex[K] fidelity[N_reviewers, K];
}

model {
  for (n in 1:N_assessments) {
    psi[n] ~ dirichlet(rep_vector(1, K));
  }
  
  for (r in 1:N_reviewers) {
    for (k in 1:K) {
      vector[K] alpha = rep_vector(1, K);
      alpha[k] = 4;
      fidelity[r, k] ~ dirichlet(alpha);
    }
  }
  
  for (n in 1:N_assessments) {
    real lds[K];
    
    // Loop over possible correct answers
    for (k in 1:K) {
      lds[k] = psi[n][k]; // Probability of kth answer being correct 
      for (r in 1:N_reviewers) {
        // Probability of response given that the kth answer is correct
        real p_response = fidelity[r, k][y[n, r]];
        lds[k] += log(p_response);
      }
    }
    target += log_sum_exp(lds);
  }
  
  for (n in 1:N_calib_assessments) {
    for (r in 1:N_reviewers) {
      target += log(fidelity[r, z_calib[n]][y_calib[n, r]]);
    }
  }
}

generated quantities {
  int<lower=0, upper=K> y_pred[N_assessments, N_reviewers];
  int<lower=0, upper=K> y_calib_pred[N_calib_assessments, N_reviewers];
    
  for (n in 1:N_assessments) {
    int z = categorical_rng(psi[n]);
    for (r in 1:N_reviewers) {
      y_pred[n, r] = categorical_rng(fidelity[r, z]);
    }
  }
  
  for (n in 1:N_calib_assessments) {
    for (r in 1:N_reviewers) {
      y_calib_pred[n, r] = categorical_rng(fidelity[r, z_calib[n]]);
    }
  }
}
