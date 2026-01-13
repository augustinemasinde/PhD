
data {
  int<lower=1> n_obs;
  array[n_obs] int<lower=0> n_pos;
  array[n_obs] int<lower=1> n_total;
  int<lower=1> age_max;
  array[n_obs] int<lower=1, upper=age_max> ages;
  matrix[n_obs, age_max] exposure_matrix;
}

parameters {
  vector[age_max] log_foi;   // No lower bound here: log can be negative
  real<lower=0, upper=1> rho;
}

transformed parameters {
  vector[n_obs] p_hat;
  for (i in 1:n_obs) {
    real sum_foi = 0;
    for (a in 1:age_max) {
      int age_diff = ages[i] - a;
      if (age_diff >= 0)
        sum_foi += exp(log_foi[a]) * pow(1 - rho, age_diff) * exposure_matrix[i, a];
    }
    p_hat[i] = 1 - exp(-sum_foi);
  }
}

model {
  // Priors
  log_foi ~ normal(-5, 1);
  rho ~ beta(1, 1);

  // Likelihood
  for (i in 1:n_obs) {
    n_pos[i] ~ binomial(n_total[i], p_hat[i]);
  }
}

generated quantities {
  vector[n_obs] log_lik;
  for (i in 1:n_obs) {
    log_lik[i] = binomial_lpmf(n_pos[i] | n_total[i], p_hat[i]);
  }
}

