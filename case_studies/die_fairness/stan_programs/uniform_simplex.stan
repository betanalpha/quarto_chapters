data {
  int<lower=1> N; // Number of rolls
}

transformed data {
  // Uniform simplex configuration
  simplex[6] upsilon = rep_vector(1.0 / 6, 6);
}

generated quantities {
  // Compute predictions
  array[6] int pred_counts = multinomial_rng(upsilon, N);
  real pred_entropy = 0;
  for (k in 1:6) {
    real q_hat = pred_counts[k] * 1.0 / N;
    pred_entropy += - lmultiply(q_hat, q_hat);
  }
}
