
data {
  int<lower=1> n_obs;
  array[n_obs] int<lower=0> n_pos;
  array[n_obs] int<lower=1> n_total;
  int<lower=1> age_max;
  array[n_obs] int<lower=0> ages;
  array[n_obs, age_max] int<lower=0, upper=1> exposure_matrix;
}
parameters {
  vector<lower=-12, upper=2>[age_max] log_foi;
  real<lower=0, upper=1> rho;
}
transformed parameters {
  vector<lower=0>[age_max] foi = exp(log_foi);
  vector[n_obs] prevalence;
  for (i in 1:n_obs) {
    real cum_hazard = 0;
    int a = ages[i];
    for (j in 1:a) {
      if (exposure_matrix[i, j] == 1) {
        cum_hazard += foi[j] * pow(1 - rho, a - j);
      }
    }
    prevalence[i] = 1 - exp(-cum_hazard);
  }
}
model {
  log_foi ~ normal(log(0.01), 1.5);
  rho ~ beta(2, 20);
  for (i in 1:n_obs) {
    real p = fmin(fmax(prevalence[i], 1e-6), 1 - 1e-6);
    n_pos[i] ~ binomial(n_total[i], p);
  }
}
generated quantities {
  vector[n_obs] log_lik;
  for (i in 1:n_obs) {
    real p = fmin(fmax(prevalence[i], 1e-6), 1 - 1e-6);
    log_lik[i] = binomial_lpmf(n_pos[i] | n_total[i], p);
  }
}

