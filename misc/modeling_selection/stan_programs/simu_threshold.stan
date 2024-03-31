data {
  int<lower=1> N; // Number of observations
}

transformed data {
  real lambda = 4.75;    // Selection threshold
  real mu = 3;           // Location of latent probability density function
  real<lower=0> tau = 2; // Shape of latent probability density function
}

generated quantities {
  // Simulated data
  real y[N];
  real N_reject = 0;
  
  for (n in 1:N) {
    // Sample latent events until one survives the selection process
    while (1) {
      y[n] = normal_rng(mu, tau);
      if (y[n] <= lambda) {
        break;
      } else {
        N_reject += 1;
      }
    }
  }
}
