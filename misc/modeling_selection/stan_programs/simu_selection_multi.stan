data {
  int<lower=1> N; // Number of multivariate observations
  int<lower=1> K; // Number of components in each observation
  
  real chi;   // Location of selection function
  real gamma; // Shape of selection function
}

generated quantities {
  // Locations of latent probability density function
  vector[K] mu;

  // Scales of latent probability density function
  vector<lower=0>[K] tau;

  // Cholesky factor of probability density function correlation matrix
  cholesky_factor_corr[K] Psi;
  
  vector[K] y[N];
  real N_reject = 0;
  
  // Simulate location and scales
  for (k in 1:K) {
    mu[k] = normal_rng(0, 5 / 2.32);
    tau[k] = fabs(normal_rng(0, 3 / 2.57));
  }
  
  // Simulate correlation matrix and compute its Cholesky decomposition
  {
    matrix[K, K] Gamma = lkj_corr_rng(K, 4 / sqrt(K));
    Psi = cholesky_decompose(Gamma);
  }
  
  for (n in 1:N) {
    // Sample latent points until one survives the selection process
    while (1) {
      vector[K] y_sample 
        = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, Psi));
      real s = sum(y_sample);
      if ( bernoulli_rng( Phi(gamma * (s - chi)) ) ) {
        y[n] = y_sample;
        break;
      } else {
        N_reject += 1;
      }
    }
  }
}
