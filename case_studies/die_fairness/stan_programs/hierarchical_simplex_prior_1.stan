transformed data {
  vector[6] ones = rep_vector(1, 6);
}

generated quantities {
  real omega = abs(normal_rng(0, 1));
  vector[6] q_baseline = dirichlet_rng(7.5 * ones);
  vector[6] q = dirichlet_rng(q_baseline / omega + ones);
}
