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
}

parameters {
  real<lower=0, upper=1> psi1; // Maximum conversion probability
  real<lower=0> psi2;          // Rate of saturation (1 / kUSD)
}

model {
  // Prior model
  psi1 ~ beta(0.5, 5.0);        // 0 <~ psi1 <~ 0.5
  psi2 ~ normal(0, 2.0 / 2.57); // 0 <~ psi2 <~ 2

  // Observational model
  for (n in 1:N) {
    target += bernoulli_lpmf(y[n] | theta(x[n], psi1, 1e-3 * psi2));
  }
}

generated quantities {
  // Posterior predictions
  array[N] real y_pred;
  real p_hat_pred;

  for (n in 1:N) {
    y_pred[n] = bernoulli_rng(theta(x[n], psi1, 1e-3 * psi2));
  }
  p_hat_pred = mean(y_pred);
}
