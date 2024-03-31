data {
  int<lower=1> N; // Number of observations
  real y[N];      // Observations
  
  int N_reject; // Number of rejected latent events
}

transformed data {
  real lambda_lower_bound = max(y);
}

parameters {
  real<lower=lambda_lower_bound> lambda; // Selection threshold
  real mu;           // Location of latent probability density function
  real<lower=0> tau; // Shape of latent probability density function
}

model {
  real log_norm = normal_lcdf(lambda | mu, tau);
  
  // Prior model
  lambda ~ normal(5, 5 / 2.32); // ~ 99% prior mass between  0 and +10
  mu ~ normal(0, 5 / 2.32);     // ~ 99% prior mass between -5 and +5
  tau ~ normal(0, 5 / 2.57);    // ~ 99% prior mass between  0 and +5
  
  // Observational model
  for (n in 1:N) {
    target += normal_lpdf(y[n] | mu, tau);
  }
  
  target += N_reject * log1m_exp(log_norm);
}

generated quantities {
  real y_pred[N];        // Posterior predictive data
  int N_reject_pred = 0; // Posterior predictive data

  for (n in 1:N) {
    // Sample latent intensities until one survives the selection process
    // Use fixed iterations instead of while loop to avoid excessive trials
    // when the selection function and latent density function are not aligned
    for (m in 1:100) {
      real y_sample = normal_rng(mu, tau);
      if (y_sample <= lambda) { 
        y_pred[n] = y_sample;
        break;
      } else {
        N_reject_pred += 1;
      }
    }
  } 
}
