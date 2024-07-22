data {
  // Observed data
  int<lower=1> N;                   // Number of unique visits
  array[N] int<lower=0, upper=1> y; // Conversion indicator
}

parameters {
  real<lower=0, upper=1> theta; // Conversion probability
}

model {
  // Prior model
  theta ~ beta(0.5, 5.0); // 0 <~ theta <~ 0.5

  // Observational model
  for (n in 1:N) {
    target += bernoulli_lpmf(y[n] | theta);
  }
}

generated quantities {
  // Posterior predictions
  array[N] real y_pred;
  real p_hat_pred;

  for (n in 1:N) {
    y_pred[n] = bernoulli_rng(theta);
  }
  p_hat_pred = mean(y_pred);
}
