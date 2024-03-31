functions {
  // Unnormalized observed probability density function
  real observed_updf(real x, real xc, real[] theta, real[] x_r, int[] x_i) {
    real chi = theta[1];
    real gamma = theta[2];
    real mu = x_r[1];
    real tau = x_r[2];
    
    return exp(normal_lpdf(x | mu, tau)) * Phi(gamma * (x - chi));
  }
}

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
  
  // Empty arrays for `integrate_1d` function
  real x_r[2] = {mu, tau};
  int x_i[0];
}

parameters {
  real chi;   // Location of selection function
  real gamma; // Shape of selection function
}

model {
  // Numerical quadrature estimator of normalization integral
  real theta[2] = {chi, gamma};
  real norm = integrate_1d(observed_updf,
                           negative_infinity(), positive_infinity(),
                           theta, x_r, x_i, 1e-8);
  
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
  real norm_error;               // Quadrature normalization errors
  real select_prob_grid[N_grid]; // Selection function values
  real y_pred[N];                // Posterior predictive data

  {
    real theta[2] = {chi, gamma};
    real quad_norm = integrate_1d(observed_updf,
                                  negative_infinity(), 
                                  positive_infinity(),
                                  theta, x_r, x_i, 1e-8);
    real exact_norm = Phi(  gamma * (mu - chi) 
                          / sqrt(1 + square(gamma * tau)) );
    norm_error = quad_norm - exact_norm;
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
