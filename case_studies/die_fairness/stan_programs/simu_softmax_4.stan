transformed data {
  matrix[6, 6] P = identity_matrix(6) + (-1.0/6.0) * rep_matrix(1, 6, 6);
  matrix[6, 5] T3 = svd_U(P)[:,1:5];
  matrix[6, 5] T4 = rep_matrix(0, 6, 5);
  matrix[6, 5] T2 = rep_matrix(0, 6, 5);

  T2[1:5,1:5] = identity_matrix(5);
  T2[6,1:5] = rep_row_vector(-1, 5);
  T4 = qr_Q(T2)[:,1:5];
}

generated quantities {
  vector[5] y;
  vector[6] x3;
  vector[6] x4;
  vector[6] q3;
  vector[6] q4;
  {
    vector[5] mu = to_vector(normal_rng(rep_array(0, 5), 0.5));
    vector[5] tau = abs(to_vector(normal_rng(rep_array(0, 5), 0.5)));
    matrix[5, 5] L = lkj_corr_cholesky_rng(5, 4.0);

    y = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, L));

    x3 = T3 * y;
    q3 = softmax(x3);

    x4 = T4 * y;
    q4 = softmax(x4);
  }
}

