data {
  int<lower=1> N; // Number of observations
  real xs[N];     // Observed x positions (m)
  real ys[N];     // Observed y positions (m)
}

transformed data {
  real L = 3; // Maximum size (m)

  // Sufficient statistics
  real min_x = min(xs);
  real max_x = max(xs);
  real min_y = min(ys);
  real max_y = max(ys);
}

parameters {
  real<lower=-L, upper=min_x> left;   // Position of left window edge (m)
  real<lower=max_x, upper=+L> right;  // Position of right window edge (m)
  real<lower=-L, upper=min_y> bottom; // Position of bottom window edge (m)
  real<lower=max_y, upper=+L> top;    // Position of top window edge (m)
}

model {
  // Implicit uniform prior density function within boundaries
  
  // Observational model
  target += - N * ( log(right - left) + log(top - bottom));
}

generated quantities {
  real min_x_pred;
  real max_x_pred;
  real min_y_pred;
  real max_y_pred;
  {
    real xs_pred[N];
    real ys_pred[N];
    for (n in 1:N) {
      xs_pred[n] = uniform_rng(left, right);
      ys_pred[n] = uniform_rng(bottom, top);
    }

    min_x_pred = min(xs_pred);
    max_x_pred = max(xs_pred);
    min_y_pred = min(ys_pred);
    max_y_pred = max(ys_pred);
  }
}
