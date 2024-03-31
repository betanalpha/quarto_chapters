functions {
  real compute_exact_norm(real chi, real gamma,
                          vector mu, vector tau, matrix Psi) {
    int K = num_elements(mu);
    real mu_s = sum(mu);
    real tau_s;

    matrix[K, K] L_Sigma = diag_pre_multiply(tau, Psi);
    matrix[K, K] Sigma = L_Sigma * L_Sigma';
    matrix[K, K] Lambda = inverse(Sigma);

    matrix[K - 1, K - 1] Gamma;
    vector[K - 1] l;

    for (i in 1:(K - 1)) {
      l[i] = Lambda[1, i + 1] - Lambda[1, 1];
      for (j in 1:(K - 1)) {
        Gamma[i, j] =  Lambda[i + 1, j + 1]
                     - Lambda[1, j + 1] - Lambda[i + 1, 1]
                     + Lambda[1, 1];
      }
    }
    tau_s = 1 / sqrt(Lambda[1, 1] - quad_form(inverse(Gamma), l));
    return Phi(  gamma * (mu_s - chi)
               / sqrt(1 + square(gamma * tau_s)) );
  }

  real compute_is_norm(vector[] y_IS, real[] ref_lpdfs,
                       real chi, real gamma,
                       vector mu, vector tau, matrix Psi) {
    int N_IS = size(y_IS);
    int K = size(y_IS[1]);

    vector[N_IS] weights;
    vector[N_IS] select_probs;
    matrix[K, K] L_Sigma = diag_pre_multiply(tau, Psi);
    for (n in 1:N_IS) {
      real active_lpdf =
        multi_normal_cholesky_lpdf(y_IS[n] | mu, L_Sigma);
      real s = sum(y_IS[n]);

      weights[n] = exp(active_lpdf - ref_lpdfs[n]);
      select_probs[n] = Phi(gamma * (s - chi));
    }
    return mean(weights .* select_probs);
  }

  // First component is the importance sampling standard error
  // Second component is the importance sampling effective sample size
  real[] compute_is_diagnostics(vector[] y_IS, real[] ref_lpdfs,
                                real chi, real gamma,
                                vector mu, vector tau, matrix Psi) {
    int N_IS = size(y_IS);
    int K = size(y_IS[1]);

    vector[N_IS] weights;
    vector[N_IS] select_probs;
    matrix[K, K] L_Sigma = diag_pre_multiply(tau, Psi);
    for (n in 1:N_IS) {
      real active_lpdf =
        multi_normal_cholesky_lpdf(y_IS[n] | mu, L_Sigma);
      weights[n] = exp(active_lpdf - ref_lpdfs[n]);
      select_probs[n] = Phi(gamma * (sum(y_IS[n]) - chi));
    }
    return { sqrt(variance(weights .* select_probs) / N_IS),
             square(sum(weights)) / dot_self(weights) };
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
  int<lower=1> N; // Number of multivariate observations
  int<lower=1> K; // Number of components in each observation
  vector[K] y[N]; // Multivariate observations

  int<lower=1> N_IS; // Size of importance sampling ensemble
}

transformed data {
  int<lower=1> M_IS = N_IS / 2; // Lower truncation for khat estimator validity

  // Generate importance sampling ensemble
  vector[K] y_IS[N_IS];
  real ref_lpdfs[N_IS];
  for (n in 1:N_IS) {
    vector[K] mu = rep_vector(0, K);
    matrix[K, K] Sigma = square(9.31) * identity_matrix(K);
    y_IS[n] = multi_normal_rng(mu, Sigma);
    ref_lpdfs[n] = multi_normal_lpdf(y_IS[n] | mu, Sigma);
  }
}

parameters {
  // Locations of latent probability density function
  vector[K] mu;

  // Scales of latent probability density function
  vector<lower=0>[K] tau;

  // Cholesky factor of probability density function correlation matrix
  cholesky_factor_corr[K] Psi;

  real chi;            // Location of selection function
  real<upper=0> gamma; // Shape of selection function
}

model {
  // Importance sampling estimator of normalization integral
  real norm = compute_is_norm(y_IS, ref_lpdfs,
                              chi, gamma, mu, tau, Psi);
  
  // Prior model
  mu ~ normal(0, 10 / 2.32);   // ~ 99% prior mass between -10 and +10
  tau ~ normal(0, 5 / 2.57);   // ~ 99% prior mass between 0 and +5
  Psi ~ lkj_corr_cholesky(4 / sqrt(K));
  chi ~ normal(0, 3 / 2.32);   // ~ 99% prior mass between -3 and +3
  gamma ~ normal(0, 3 / 2.57); // ~ 99% prior mass between -3 and  0
  
  // Observational model
  for (n in 1:N) {
    real s = sum(y[n]);
    real log_select_prob = log( Phi(gamma * (s - chi)) );
    real lpdf =
      multi_normal_cholesky_lpdf(y[n] | mu, diag_pre_multiply(tau, Psi));
    target += log_select_prob + lpdf - log(norm);
  }
}

generated quantities {
  real norm_error; // Importance sampling normalization errors
  real ISSE;       // Importance sampling standard error
  real ISESS;      // Importance sampling effective sample size

  vector[N_IS] weights; // Importance sampling weights
  real khat;            // khat statistic of importance sampling weights

  vector[K] y_pred[N]; // Posterior predictive data
  
  {
    real exact_norm = compute_exact_norm(chi, gamma, mu, tau, Psi);
    real norm = compute_is_norm(y_IS, ref_lpdfs,
                                chi, gamma, mu, tau, Psi);
    real diagnostics[2] = compute_is_diagnostics(y_IS, ref_lpdfs,
                                                 chi, gamma, mu, tau, Psi);

    norm_error = norm - exact_norm;
    ISSE = diagnostics[1];
    ISESS = diagnostics[2];

    for (n in 1:N_IS) {
      real active_lpdf =
        multi_normal_cholesky_lpdf(y_IS[n] | mu, diag_pre_multiply(tau, Psi));
      weights[n] = exp(active_lpdf - ref_lpdfs[n]);
    }
    khat = compute_khat(weights[M_IS:N_IS]);
  }

  for (n in 1:N) {
    // Sample latent intensities until one survives the selection
    // process.  Use fixed iterations instead of while loop to avoid 
    // excessive trials when the selection function and latent density 
    // function are not aligned.
    for (m in 1:100) {
      vector[K] y_sample
        = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, Psi));
      real s = sum(y_sample);
      if ( bernoulli_rng( Phi(gamma * (s - chi)) ) ) {
        y_pred[n] = y_sample;
        break;
      }
    }
  }
}
