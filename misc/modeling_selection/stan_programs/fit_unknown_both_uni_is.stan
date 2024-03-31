functions {
real compute_is_norm(real[] y_IS, real[] ref_lpdfs,
                       real chi, real gamma, real mu, real tau) {
    int N_IS = size(y_IS);
    vector[N_IS] weights;
    vector[N_IS] select_probs;
    for (n in 1:N_IS) {
      weights[n] = exp(normal_lpdf(y_IS[n] | mu, tau) - ref_lpdfs[n]);
      select_probs[n] = Phi(gamma * (y_IS[n] - chi));
    }
    return mean(weights .* select_probs);
  }

  // First returned component is the importance sampling standard error
  // Second returned component is the importance sampling effective sample size
  real[] compute_is_diagnostics(real[] y_IS, real[] ref_lpdfs,
                                real chi, real gamma,
                                real mu, real tau) {
    int N_IS = size(y_IS);
    vector[N_IS] weights;
    vector[N_IS] select_probs;
    for (n in 1:N_IS) {
      weights[n] = exp(normal_lpdf(y_IS[n] | mu, tau) - ref_lpdfs[n]);
      select_probs[n] = Phi(gamma * (y_IS[n] - chi));
    }
    return {sqrt(variance(weights .* select_probs) / N_IS),
            square(sum(weights)) / dot_self(weights)};
  }

  real compute_khat(vector fs) {
    int N = num_elements(fs);
    vector[N] sorted_fs = sort_asc(fs);

    real R;
    real q;
    real M;
    vector[N] b_hat_vec;
    vector[N] log_w_vec;

    real khat;
    real max_log_w;
    real b_hat_denom = 0;
    real b_hat_numer = 0;
    real b_hat;

    if (sorted_fs[1] == sorted_fs[N]) return not_a_number();
    if (sorted_fs[1] < 0) return not_a_number();

    // Estimate 25% quantile
    R = floor(0.25 * N + 0.5);
    for (n in 1:N) {
      if (n + 0.5 >= R) {
        q = sorted_fs[n];
        break;
      }
    }
    if (q == sorted_fs[1]) return not_a_number();

    // Heuristic Pareto configuration
    M = 20 + floor(sqrt(N));

    for (m in 1:N) {
      if (m > M) break;

      b_hat_vec[m] = 1 / sorted_fs[N] + (1 - sqrt(M / (m - 0.5))) / (3 * q);
      if (b_hat_vec[m] != 0) {
        khat = - mean( log(1 - (b_hat_vec[m] * sorted_fs) ) );
        log_w_vec[m] = N * ( log(b_hat_vec[m] / khat) + khat - 1);
      } else {
        log_w_vec[m] = 0;
      }
    }

    max_log_w = log_w_vec[1];
    for (m in 1:N) {
      if (m > M) break;
      if (log_w_vec[m] > max_log_w) max_log_w = log_w_vec[m];
    }

    for (m in 1:N) {
      if (m <= M) {
        real weight = exp(log_w_vec[m] - max_log_w);
        b_hat_numer += b_hat_vec[m] * weight;
        b_hat_denom += weight;
      } else {
        break;
      }
    }
    b_hat = b_hat_numer / b_hat_denom;

    return -mean( log(1 - b_hat * sorted_fs) );
  }
}

data {
  int<lower=1> N; // Number of observations
  real y[N];      // Observations

  int<lower=1> N_grid; // Number of latent value grid points
  real y_grid[N_grid]; // Latent value grid points

  int<lower=1> N_IS; // Size of importance sampling ensemble
}

transformed data {
  int<lower=1> M_IS = N_IS / 2; // Lower truncation for khat estimator validity

  // Generate importance sampling ensemble
  real y_IS[N_IS];
  real ref_lpdfs[N_IS];
  for (n in 1:N_IS) {
    y_IS[n] = normal_rng(0, 7.2);
    ref_lpdfs[n] = normal_lpdf(y_IS[n] | 0, 7.2);
  }
}

parameters {
  real mu;           // Location of latent probability density function
  real<lower=0> tau; // Shape of latent probability density function
  real chi;   // Location of selection function
  real gamma; // Shape of selection function
}

model {
  // Importance sampling estimator of normalization integral
  real norm = compute_is_norm(y_IS, ref_lpdfs, chi, gamma, mu, tau);
  
  // Prior model
  mu ~ normal(0, 5 / 2.32);   // ~ 99% prior mass between -5 and +5
  tau ~ normal(0, 5 / 2.57);  // ~ 99% prior mass between 0 and +5
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
  real norm_error; // Importance sampling normalization errors
  real ISSE;       // Importance sampling standard error
  real ISESS;      // Importance sampling effective sample size

  vector[N_IS] weights; // Importance sampling weights
  real khat;            // khat statistic of importance sampling weights
  
  real select_prob_grid[N_grid]; // Selection function values
  real y_pred[N];                // Posterior predictive data
  
{
    real exact_norm = Phi(  gamma * (mu - chi)
                          / sqrt(1 + square(gamma * tau)) );
    real norm = compute_is_norm(y_IS, ref_lpdfs,
                                chi, gamma, mu, tau);
    real diagnostics[2] = compute_is_diagnostics(y_IS, ref_lpdfs,
                                                 chi, gamma, mu, tau);

    norm_error = norm - exact_norm;
    ISSE = diagnostics[1];
    ISESS = diagnostics[2];

    for (n in 1:N_IS)
      weights[n] = exp(normal_lpdf(y_IS[n] | mu, tau) - ref_lpdfs[n]);
    khat = compute_khat(weights[M_IS:N_IS]);
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
