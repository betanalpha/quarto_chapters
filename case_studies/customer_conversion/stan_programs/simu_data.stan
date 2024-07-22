functions {
  // Saturating conversion probability function
  real theta(real x, real psi1, real psi2) {
    return psi1 * (-expm1(-psi2 * x));
  }
}

transformed data {
  real<lower=0, upper=1> psi1 = 0.47;
  real<lower=0> psi2 = 0.5 / 1000;
  real<lower=0, upper=1> lambda = 0.3;
  real<lower=0, upper=1> theta_VIP = 0.95;

  int<lower=1> N = 1500;
  int<lower=1> N_aux = 200;
}

generated quantities {
  array[N] int<lower=0, upper=1> y;
  array[N] real<lower=0> x;

  array[N_aux] int<lower=0, upper=1> y_aux;

  for (n in 1:N) {
    x[n] = fabs(normal_rng(0, 4500 / 2.57));
    if (bernoulli_rng(lambda)) {
      y[n] = bernoulli_rng(theta_VIP);
    } else {
      y[n] = bernoulli_rng(theta(x[n], psi1, psi2));
    }
  }

  for (n in 1:N_aux) {
    y_aux[n] = bernoulli_rng(theta_VIP);
  }
}
