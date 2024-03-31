data {
  int<lower=1> N; // Number of multivariate observations
  int<lower=1> K; // Number of components in each observation
  vector[K] y[N]; // Multivariate observations
}

parameters {
  // Locations of latent probability density function
  vector[K] mu;

  // Scales of latent probability density function
  vector<lower=0>[K] tau;

  // Cholesky factor of probability density function correlation matrix
  cholesky_factor_corr[K] Psi;
}

model {
  // Prior model
  mu ~ normal(0, 10 / 2.32);  // ~ 99% prior mass between -10 and +10
  tau ~ normal(0, 5 / 2.57);  // ~ 99% prior mass between 0 and +5
  Psi ~ lkj_corr_cholesky(4 / sqrt(K));
  
  // Observational model
  for (n in 1:N) {
    target += multi_normal_cholesky_lpdf(y[n] | mu,
                                                diag_pre_multiply(tau, Psi));
  }
}

generated quantities {
  vector[K] y_pred[N]; // Posterior predictive data
  
  for (n in 1:N) {
    y_pred[n] = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, Psi));
  }
}
