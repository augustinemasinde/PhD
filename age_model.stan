
data {
  int<lower=1> N;
  array[N] int<lower=1> ages;
  array[N] int<lower=0> positive;
  array[N] int<lower=0> total;
}

parameters {
  vector<lower=0>[N] foi;  // force of infection per age group, positive
  real<lower=0, upper=1> rho; // reporting probability
}

model {
  // Priors
  foi ~ normal(0, 1);
  rho ~ beta(2, 2);

  // Likelihood: positive ~ binomial(total, p_hat)
  for (n in 1:N) {
    real p_hat = 1 - exp(-foi[n] * ages[n]) * (1 - rho);
    positive[n] ~ binomial(total[n], p_hat);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real p_hat = 1 - exp(-foi[n] * ages[n]) * (1 - rho);
    log_lik[n] = binomial_lpmf(positive[n] | total[n], p_hat);
  }
}

