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

  // Number of entrants in each race who did not finish
  array[N_races] int<lower=0, upper=N_entrants> race_N_entrants_dnf;

  // Indices for extracting did not finish entrant information in each race
  array[N_races] int race_dnf_start_idxs;
  array[N_races] int race_dnf_end_idxs;

  // Total number of finishes across all races
  int <lower=1> N_entrances_fs;
  array[N_entrances_fs] int race_entrant_f_idxs;   // Entrant index
  array[N_entrances_fs] real race_entrant_f_times; // Entrant finish times

  // Total number of forfeits across all races
  int<lower=0> N_entrances_dnfs;
  array[N_entrances_dnfs] int race_entrant_dnf_idxs; // Entrant index

  // MapRando versioning
  int<lower=1> N_versions;
  array[N_races] int<lower=1, upper=N_versions> version_idxs;
}

parameters {
  real gamma;                       // Log baseline finish time (log seconds)

  vector[N_races] eta_difficulties; // Non-centered seed difficulties
  vector<lower=0>[N_versions] tau_difficulties; // Seed difficulty population scale

  vector[N_entrants] eta_skills;    // Non-centered entrant skills
  real<lower=0> tau_skills;         // Entrant skill population scale

  real<lower=0> psi;                // Inverse gamma dispersion configuration

  array[N_entrants] real kappas;         // Forfeit thresholds
  array[N_entrants] real<lower=0> betas; // Forfeit scalings
}

transformed parameters {
  // Derive centered race difficulties and entrant skills
  vector[N_races] difficulties
    = tau_difficulties[version_idxs] .* eta_difficulties;
  vector[N_entrants] skills
    = tau_skills * eta_skills;
}

model {
  // Prior model
  gamma ~ normal(8.045, 0.237); // log(1800 s) < gamma < log(5400 s)
  eta_difficulties ~ normal(0, 1);     // Non-centered individual model
  tau_difficulties ~ normal(0, 0.270); // 0 <~ tau_difficulties <~ log(2)
  eta_skills ~ normal(0, 1);           // Non-centered individual model
  tau_skills ~ normal(0, 0.270);       // 0 <~    tau_skills    <~ log(2)
  psi ~ normal(0, 0.389);       // 0 <~ psi <~ 1
  kappas ~ normal(0, 2.16);     // -5 <~ kappas <~ +5
  betas ~ normal(0, 1.95);      // 0 <~ betas <~ +5

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
      real delta = difficulties[r] - skills[entrant_idx];
      real mu = exp(gamma + delta);
      real logit_q = betas[entrant_idx] * (delta - kappas[entrant_idx]);

      entrant_f_times[n] ~ inv_gamma_md(mu, psi);
      0 ~ bernoulli_logit(logit_q); // Did not forfeit
    }

    if (race_N_entrants_dnf[r] > 0) {
      // Extract details for entrants who forfeited
      int N_entrants_dnf = race_N_entrants_dnf[r];
      array[N_entrants_dnf]
        int dnf_idxs = linspaced_int_array(N_entrants_dnf,
                                           race_dnf_start_idxs[r],
                                           race_dnf_end_idxs[r]);
      array[N_entrants_dnf]
        int entrant_dnf_idxs = race_entrant_dnf_idxs[dnf_idxs];

      // Forfeited entrant model
      for (n in 1:N_entrants_dnf) {
        int entrant_idx = entrant_dnf_idxs[n];
        real delta = difficulties[r] - skills[entrant_idx];
        real logit_q = betas[entrant_idx] * (delta - kappas[entrant_idx]);

        1 ~ bernoulli_logit(logit_q); // Did forfeit
      }
    }
  }
}

generated quantities {
  // Posterior predictions
  array[N_entrances_fs] real race_entrant_f_times_pred;
  array[N_races] int race_N_entrants_dnf_pred = rep_array(0, N_races);

  for (r in 1:N_races) {
    // Extract details for entrants who finished
    int N_entrants_f = race_N_entrants_f[r];
    array[N_entrants_f] int f_idxs = linspaced_int_array(N_entrants_f,
                                                         race_f_start_idxs[r],
                                                         race_f_end_idxs[r]);
    array[N_entrants_f] int entrant_f_idxs = race_entrant_f_idxs[f_idxs];

    for (n in 1:N_entrants_f) {
      int entrant_idx = entrant_f_idxs[n];
      real delta = difficulties[r] - skills[entrant_idx];
      real mu = exp(gamma + delta);
      real logit_q = betas[entrant_idx] * (delta - kappas[entrant_idx]);

      // Finish time prediction conditioned on not forfeiting
      race_entrant_f_times_pred[race_f_start_idxs[r] + n - 1]
        = inv_gamma_md_rng(mu, psi);

      // Forfeit prediction
      race_N_entrants_dnf_pred[r] += bernoulli_logit_rng(logit_q);
    }

    if (race_N_entrants_dnf[r] > 0) {
      // Extract details for entrants who forfeited
      int N_entrants_dnf = race_N_entrants_dnf[r];
      array[N_entrants_dnf]
        int dnf_idxs = linspaced_int_array(N_entrants_dnf,
                                           race_dnf_start_idxs[r],
                                           race_dnf_end_idxs[r]);
      array[N_entrants_dnf]
        int entrant_dnf_idxs = race_entrant_dnf_idxs[dnf_idxs];

      for (n in 1:N_entrants_dnf) {
        int entrant_idx = entrant_dnf_idxs[n];
        real delta = difficulties[r] - skills[entrant_idx];
        real logit_q = betas[entrant_idx] * (delta - kappas[entrant_idx]);

        // Forfeit prediction
        race_N_entrants_dnf_pred[r] += bernoulli_logit_rng(logit_q);
      }
    }
  }
}
