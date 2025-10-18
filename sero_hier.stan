
data {
  int<lower=1> I;
  int<lower=0> Npos[I];
  int<lower=0> Ntest[I];
  vector[I] w;
  int<lower=1> S;
  int<lower=1> A;
  int<lower=1> G;
  int<lower=1,upper=S> site[I];
  int<lower=1,upper=A> age[I];
  int<lower=1,upper=G> sex[I];
}
parameters {
  real mu;
  vector[S] alpha_site;
  vector[A] alpha_age;
  vector[G] alpha_sex;
  real<lower=0,upper=1> sens;
  real<lower=0,upper=1> spec;
}
transformed parameters {
  vector[I] p_true;
  vector[I] p_obs;
  for (i in 1:I) {
    real lin = mu + alpha_site[site[i]] + alpha_age[age[i]] + alpha_sex[sex[i]];
    p_true[i] = inv_logit(lin);
    p_obs[i]  = p_true[i] * sens + (1 - p_true[i]) * (1 - spec);
  }
}
model {
  // Weakly informative priors
  mu ~ normal(0, 1.5);
  alpha_site ~ normal(0, 1);
  alpha_age  ~ normal(0, 1);
  alpha_sex  ~ normal(0, 1);

  // Priors for test characteristics: replace with your preferred Beta(a,b)
  spec ~ beta(2765.505, 194.95);
  sens ~ beta(1321.92,55.08);

  // Likelihood
  for (i in 1:I)
    Npos[i] ~ binomial(Ntest[i], p_obs[i]);
}
generated quantities {
  real pop_prev = 0;
  for (i in 1:I) pop_prev = pop_prev + w[i] * p_true[i];
  vector[I] p_stratum = p_true;
}

