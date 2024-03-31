data {
  int<lower=1> N; // Number of observations
}

transformed data {
  real mu = -1;          // Location of latent probability density function
  real<lower=0> tau = 3; // Shape of latent probability density function
  real chi = 2;          // Location of selection function
  real gamma = 0.75;     // Shape of selection function
}

generated quantities {
  // Simulate filtered data
  real y[N];
  real N_reject = 0;
  
  for (n in 1:N) {
    // Sample latent points until one survives the selection process
    while (1) {
      y[n] = normal_rng(mu, tau);
      if ( bernoulli_rng( Phi(gamma * (y[n] - chi)) ) ) {
        break;
      } else {
        N_reject += 1;
      }
    }
  }
}
