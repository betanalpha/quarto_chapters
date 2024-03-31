functions {
  real compute_norm(real chi, real gamma, 
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
}

data {
  int<lower=1> N; // Number of multivariate observations
  int<lower=1> K; // Number of components in each observation
  vector[K] y[N]; // Multivariate observations
  
  // Locations of latent probability density function
  vector[K] mu;

  // Scales of latent probability density function
  vector<lower=0>[K] tau;

  // Cholesky factor of probability density function correlation matrix
  cholesky_factor_corr[K] Psi;
}

parameters {
  real chi;            // Location of selection function
  real<upper=0> gamma; // Shape of selection function
}

model {
  real norm = compute_norm(chi, gamma, mu, tau, Psi);
  
  // Prior model
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
  vector[K] y_pred[N]; // Posterior predictive data
  
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
