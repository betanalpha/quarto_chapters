data {
  int<lower=1> N;                         // Number of rolls
  array[N] int<lower=1, upper=6> outcome; // Roll outcomes

  int<lower=1> N_dice;                          // Number of dice
  array[N] int<lower=1, upper=N_dice> die_idxs; // Die rolled
}

transformed data {
  // Compute outcome totals
  array[6] int counts = rep_array(0, 6);
  array[N_dice, 6] int counts_per_die;
  array[N_dice] int N_per_die;
  
  for (n in 1:N) {
    counts[outcome[n]] += 1;
  }

  for (i in 1:N_dice) {
    counts_per_die[i] = rep_array(0, 6);
    for (n in 1:N) {
      counts_per_die[die_idxs[n], outcome[n]] += 1;
    }
    N_per_die[i] = sum(counts_per_die[i]);
  }
}

parameters {
  array[N_dice] simplex[6] q;
}

model {
  for (i in 1:N_dice) {
    q[i] ~ dirichlet(rep_vector(1, 6));
    counts_per_die[i] ~ multinomial(q[i]);
  }
}

generated quantities {
  array[6] int pred_counts = rep_array(0, 6);
  real pred_entropy = 0;

  array[N_dice, 6] real pred_freq_die;
  array[N_dice] real pred_entropy_die;

  array[6] real pred_var_freq;
  real pred_var_entropy;

  for (i in 1:N_dice) {
    array[6] int pred_counts_die = multinomial_rng(q[i], N_per_die[i]);

    pred_entropy_die[i] = 0;
    for (k in 1:6) {
      pred_counts[k] += pred_counts_die[k];
      pred_freq_die[i, k] = pred_counts_die[k] * 1.0 / N_per_die[i];
      pred_entropy_die[i] += - lmultiply(pred_freq_die[i, k],
                                         pred_freq_die[i, k]);
    }
  }

  for (k in 1:6) {
    real q_hat = pred_counts[k] * 1.0 / N;
    pred_entropy += - lmultiply(q_hat, q_hat);

    pred_var_freq[k] = variance(pred_freq_die[:,k]);
  }

  pred_var_entropy = variance(pred_entropy_die);
}
