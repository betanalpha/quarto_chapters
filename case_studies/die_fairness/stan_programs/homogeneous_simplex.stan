data {
  int<lower=1> N;                         // Number of rolls
  array[N] int<lower=1, upper=6> outcome; // Roll outcomes
}

transformed data {
  // Compute outcome totals
  array[6] int counts = rep_array(0, 6);
  for (n in 1:N) {
    counts[outcome[n]] += 1;
  }
}

parameters {
  // General simplex configuration
  simplex[6] q;
}

model {
  q ~ dirichlet(rep_vector(1, 6));
  counts ~ multinomial(q);
}

generated quantities {
  // Posterior predictions
  array[6] int pred_counts = multinomial_rng(q, N);
  real pred_entropy = 0;

  for (k in 1:6) {
    real q_hat = pred_counts[k] * 1.0 / N;
    pred_entropy += - lmultiply(q_hat, q_hat);
  }
}
