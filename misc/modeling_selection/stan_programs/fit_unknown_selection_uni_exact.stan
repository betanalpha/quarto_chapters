data {
  int<lower=1> N; // Number of observations
  real y[N];      // Observations

  int<lower=1> N_grid; // Number of latent value grid points
  real y_grid[N_grid]; // Latent value grid points
}

transformed data {
  // Exact latent probability distribution configuration
  real mu = -1;          // Location
  real<lower=0> tau = 3; // Shape
}

parameters {
  real chi;   // Location of selection function
  real gamma; // Shape of selection function
}

model {
  // Filtered normalization
  real norm = Phi( gamma * (mu - chi) / sqrt(1 + square(gamma * tau)) );
  
  // Prior model
  chi ~ normal(0, 3 / 2.32);   // ~ 99% prior mass between -3 and +3
  gamma ~ normal(0, 3 / 2.32); // ~ 99% prior mass between -3 and +3
  
  // Observational model
  for (n in 1:N) {
    real log_select_prob = log( Phi(gamma * (y[n] - chi)) );
    real lpdf = normal_lpdf(y[n] | mu, tau);
    target += log_select_prob + lpdf - log(norm);
  }
}

generated quantities {
  real select_prob_grid[N_grid]; // Selection function values
  real y_pred[N];                // Posterior predictive data
  
  for (n in 1:N_grid)
    select_prob_grid[n] = Phi(gamma * (y_grid[n] - chi));
  
  for (n in 1:N) {
    // Sample latent intensities until one survives the selection
    // process.  Use fixed iterations instead of while loop to avoid 
    // excessive trials when the selection function and latent density 
    // function are not aligned.
    for (m in 1:100) {
      real y_sample = normal_rng(mu, tau);
      if ( bernoulli_rng( Phi(gamma * (y_sample - chi)) ) ) {
        y_pred[n] = y_sample;
        break;
      }
    }
  }
}
