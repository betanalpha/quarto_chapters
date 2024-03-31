data {
  int<lower=1> N; // Number of observations
  real y[N];      // Observations
}

parameters {
  real mu;           // Location of latent probability density function
  real<lower=0> tau; // Shape of latent probability density function
}

model {
  // Prior model
  mu ~ normal(0, 5 / 2.32);     // ~ 99% prior mass between -5 and +5
  tau ~ normal(0, 5 / 2.57);    // ~ 99% prior mass between  0 and +5
  
  // Observational model
  target += normal_lpdf(y | mu, tau);
}

generated quantities {
  real y_pred[N]; // Posterior predictive data

  for (n in 1:N) {
    y_pred[n] = normal_rng(mu, tau);
  } 
}
