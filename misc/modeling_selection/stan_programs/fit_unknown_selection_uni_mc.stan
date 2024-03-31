functions {
  real compute_mc_norm(real[] y_MC, real chi, real gamma) {
    int N_MC = size(y_MC);
    real select_probs[N_MC];
    for (n in 1:N_MC) 
      select_probs[n] = Phi(gamma * (y_MC[n] - chi));
    return mean(select_probs);
  }
  real compute_MCSE(real[] y_MC, real chi, real gamma) {
    int N_MC = size(y_MC);
    real select_probs[N_MC];
    for (n in 1:N_MC) 
      select_probs[n] = Phi(gamma * (y_MC[n] - chi));
    return sqrt(variance(select_probs) / N_MC);
  }
}

data {
  int<lower=1> N; // Number of observations
  real y[N];      // Observations

  int<lower=1> N_grid; // Number of latent value grid points
  real y_grid[N_grid]; // Latent value grid points

  int<lower=1> N_MC; // Size of Monte Carlo ensemble
}

transformed data {
  // Exact latent probability distribution configuration
  real mu = -1;          // Location
  real<lower=0> tau = 3; // Shape
  
  // Generate Monte Carlo ensemble
  real y_MC[N_MC];
  for (n in 1:N_MC) 
    y_MC[n] = normal_rng(mu, tau);
}

parameters {
  real chi;   // Location of selection function
  real gamma; // Shape of selection function
}

model {
  // Monte Carlo estimator of normalization integral
  real norm = compute_mc_norm(y_MC, chi, gamma);
  
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
  // Exact Monte Carlo error
  real norm_error;
  // Monte Carlo standard errors for normalization estimation
  real MCSE = compute_MCSE(y_MC, chi, gamma);
  
  real select_prob_grid[N_grid]; // Selection function values
  real y_pred[N];                // Posterior predictive data
  
  {
    real exact_norm = Phi(  gamma * (mu - chi) 
                          / sqrt(1 + square(gamma * tau)) );
    real norm = compute_mc_norm(y_MC, chi, gamma);
    norm_error = norm - exact_norm;
  }
  
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
