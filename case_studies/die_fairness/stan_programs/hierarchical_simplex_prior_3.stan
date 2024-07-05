transformed data {
  matrix[6, 6] P = identity_matrix(6) + (-1.0/6.0) * rep_matrix(1, 6, 6);
  matrix[6, 5] T3 = svd_U(P)[:,1:5];
}

generated quantities {
  vector[6] q_baseline;
  vector[6] q;
  {
    vector[5] mu = to_vector(normal_rng(rep_array(0, 5), 0.33));
    vector[5] tau = abs(to_vector(normal_rng(rep_array(0, 5), 0.66)));
    matrix[5, 5] L = lkj_corr_cholesky_rng(5, 4.0);

    vector[5] y = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, L));

    q_baseline = softmax(T3 * mu);
    q = softmax(T3 * y);
  }
}
