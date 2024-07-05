transformed data {
  matrix[6, 5] T1 = rep_matrix(0, 6, 5);
  T1[1:5,1:5] = identity_matrix(5);
}

generated quantities {
  vector[5] y;
  vector[6] x;
  vector[6] q;
  {
    vector[5] mu = to_vector(normal_rng(rep_array(0, 5), 0.5));
    vector[5] tau = abs(to_vector(normal_rng(rep_array(0, 5), 0.5)));
    matrix[5, 5] L = lkj_corr_cholesky_rng(5, 4.0);

    y = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, L));
    x = T1 * y;
    q = softmax(x);
  }
}
