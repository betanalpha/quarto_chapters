generated quantities {
  vector[6] x;
  vector[6] q;
  {
    vector[6] mu = to_vector(normal_rng(rep_array(0, 6), 0.5));
    vector[6] tau = abs(to_vector(normal_rng(rep_array(0, 6), 0.5)));
    matrix[6, 6] L = lkj_corr_cholesky_rng(6, 4.0);

    x = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, L));
    q = softmax(x);
  }
}
