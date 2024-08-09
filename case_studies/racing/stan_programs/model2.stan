functions {
  // Mean-dispersion parameterization of inverse gamma family
  real inv_gamma_md_lpdf(real x, real mu, real psi) {
    return inv_gamma_lpdf(x | inv(psi) + 2, mu * (inv(psi) + 1));
  }

  real inv_gamma_md_rng(real mu, real psi) {
    return inv_gamma_rng(inv(psi) + 2, mu * (inv(psi) + 1));
  }
}

data {
  int<lower=1> N_races;    // Total number of races
  int<lower=1> N_entrants; // Total number of entrants
  // Each entrant is assigned a unique index in [1, N_entrants]

  // Number of entrants in each race who finished
  array[N_races] int<lower=1, upper=N_entrants> race_N_entrants_f;

  // Indices for extracting finished entrant information in each race
  array[N_races] int race_f_start_idxs;
  array[N_races] int race_f_end_idxs;

  // Total number of entrant finishes across all races
  int <lower=1> N_entrances_fs;

  // Finished entrant indices within each race
  array[N_entrances_fs] int race_entrant_f_idxs;

  // Entrant finish times within each race
  array[N_entrances_fs] real race_entrant_f_times;
}

parameters {
  real gamma;                       // Log baseline finish time (log seconds)
  array[N_races] real difficulties; // Seed difficulties
  array[N_entrants] real skills;    // Entrant skills
  real<lower=0> psi;                // Inverse gamma dispersion configuration
}

model {
  // Prior model
  gamma ~ normal(8.045, 0.237);    // log(1800 s) < gamma < log(5400 s)
  difficulties ~ normal(0, 0.299); // -log(2) <~ difficulties <~ +log(2)
  skills ~ normal(0, 0.299);       // -log(2) <~    skills    <~ +log(2)
  psi ~ normal(0, 0.389);          // 0 <~ psi <~ 1

  // Observational model
  for (r in 1:N_races) {
    // Extract details for entrants who finished
    int N_entrants_f = race_N_entrants_f[r];
    array[N_entrants_f] int f_idxs = linspaced_int_array(N_entrants_f,
                                                         race_f_start_idxs[r],
                                                         race_f_end_idxs[r]);
    array[N_entrants_f] int entrant_f_idxs = race_entrant_f_idxs[f_idxs];
    array[N_entrants_f] real entrant_f_times = race_entrant_f_times[f_idxs];

    // Finished entrant model
    for (n in 1:N_entrants_f) {
      int entrant_idx = entrant_f_idxs[n];
      real mu = exp(gamma + difficulties[r] - skills[entrant_idx]);
      entrant_f_times[n] ~ inv_gamma_md(mu, psi);
    }
  }
}

generated quantities {
  // Posterior predictions
  array[N_entrances_fs] real race_entrant_f_times_pred;

  for (r in 1:N_races) {
    // Extract details for entrants who finished
    int N_entrants_f = race_N_entrants_f[r];
    array[N_entrants_f] int f_idxs = linspaced_int_array(N_entrants_f,
                                                         race_f_start_idxs[r],
                                                         race_f_end_idxs[r]);
    array[N_entrants_f] int entrant_f_idxs = race_entrant_f_idxs[f_idxs];

    // Finish time predictions conditioned on not forfeiting
    for (n in 1:N_entrants_f) {
      int entrant_idx = entrant_f_idxs[n];
      real mu = exp(gamma + difficulties[r] - skills[entrant_idx]);
      race_entrant_f_times_pred[race_f_start_idxs[r] + n - 1]
        = inv_gamma_md_rng(mu, psi);
    }
  }
}
