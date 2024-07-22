functions {
  // Saturating conversion probability function
  real theta(real x, real psi1, real psi2) {
    return psi1 * (-expm1(-psi2 * x));
  }
}

data {
  // Observed data
  int<lower=1> N;                    // Number of unique visits
  array[N] int<lower=0, upper=1> y;  // Conversion indicator
  array[N] real x;                   // Previous purchases (USD)

  int<lower=1> N_aux;                       // Number of VIP-only visits
  array[N_aux] int<lower=0, upper=1> y_aux; // Conversion indicator
}

parameters {
  real<lower=0, upper=1> psi1;      // Maximum conversion probability
  real<lower=0> psi2;               // Rate of saturation (1 / kUSD)
  real<lower=0, upper=1> lambda;    // Proportion of VIP visitors
  real<lower=0, upper=1> theta_VIP; // VIP conversion probability
}

model {
  // Prior model
  psi1 ~ beta(0.5, 5.0);        // 0   <~    psi1   <~ 0.5
  psi2 ~ normal(0, 2.0 / 2.57); // 0   <~    psi2   <~ 2
  lambda ~ beta(1, 1);          // Uniform probability density function
  theta_VIP ~ beta(5.0, 0.5);   // 0.5 <~ theta_VIP <~ 1

  // Observational model
  for (n in 1:N) {
    target += log_mix(lambda,
                      bernoulli_lpmf(y[n] | theta_VIP),
                      bernoulli_lpmf(y[n] | theta(x[n], psi1, 1e-3 * psi2)));
  }

  for (n in 1:N_aux) {
    target += bernoulli_lpmf(y_aux[n] | theta_VIP);
  }
}

generated quantities {
  // Posterior predictions
  array[N] real y_pred;
  real p_hat_pred;

  array[N_aux] real y_aux_pred;
  real p_hat_aux_pred;

  for (n in 1:N) {
    if (bernoulli_rng(lambda)) {
      y_pred[n] = bernoulli_rng(theta_VIP);
    } else {
      y_pred[n] = bernoulli_rng(theta(x[n], psi1, 1e-3 * psi2));
    }
  }
  p_hat_pred = mean(y_pred);

  for (n in 1:N_aux) {
    y_aux_pred[n] = bernoulli_rng(theta_VIP);
  }
  p_hat_aux_pred = mean(y_aux_pred);
}
