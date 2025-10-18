data {
  int <lower=0> A;                     // number of age classes
  int <lower=0> A2;                    // number of FOI groups
  int <lower=0> N;                     // number of individuals
  int <lower=0> age[N];                // age for each individual
  int <lower=0, upper=1> Y[N];         // outcome (0/1) for each individual
  int <lower=1, upper=A2> ind_by_age[A]; // FOI group index by age

  int <lower=0> NTESH;                 // number of individuals in Tesh
  int <lower=0> NGTESH;                // number of age groups in Tesh
  int <lower=0> SizeGTESH[NGTESH];     // size (years) of each age group
  int <lower=0> GTESH[NTESH];          // age group of each Tesh individual
  int <lower=0, upper=1> YTESH[NTESH]; // outcome (0/1) for Tesh individuals
  int <lower=0> GageTESH[A];           // FOI group of each age
  int <lower=0> Y1;                    // year 1 (start of fixed FOI period)
  int <lower=0> YearsFixed;            // number of years with fixed FOI
}

parameters {
  real logitlambdaGP[A2];              // logit FOI by FOI group
}

transformed parameters {
  real lambda[A];                      // FOI by age
  real cum_foi[A];                     // cumulative FOI by age
  real cum_foi2[A];                    // cumulative FOI for Tesh lot
  real logitlambda[A];                 // logit FOI by age
  real logitlambda2[A];                // logit FOI adjusted for fixed period
  real meancumfoiG[NGTESH];            // mean cumulative FOI by Tesh group
  real sumcumfoiG[NGTESH];             // sum FOI by Tesh group

  // enforce minimum logit value
  for (j in 1:A) {
    if (logitlambdaGP[ind_by_age[j]] < -6)
      logitlambda[j] = -6;
    else
      logitlambda[j] = logitlambdaGP[ind_by_age[j]];
  }

  // fix lambda for recent years
  for (k in 1:A) {
    if ((A - YearsFixed) < k)
      logitlambda2[k] = logitlambda[A - YearsFixed];
    else
      logitlambda2[k] = logitlambda[k];
  }

  // transform to probability scale
  for (j in 1:A)
    lambda[j] = inv_logit(logitlambda2[j]);

  // cumulative FOI
  cum_foi[1] = lambda[1];
  for (j in 2:A)
    cum_foi[j] = cum_foi[j-1] + lambda[j];

  // cumulative FOI adjusted for start year
  for (j in 1:A) {
    if (j < Y1)
      cum_foi2[j] = 0;
    else
      cum_foi2[j] = cum_foi[j] - cum_foi[Y1];
  }

  // initialize sums
  for (j in 1:NGTESH)
    sumcumfoiG[j] = 0;

  // sum FOI by Tesh group
  for (j in 1:A) {
    if (GageTESH[j] != 0)
      sumcumfoiG[GageTESH[j]] += cum_foi2[j];
  }

  // average FOI by Tesh group
  for (j in 1:NGTESH)
    meancumfoiG[j] = sumcumfoiG[j] / SizeGTESH[j];
}

model {
  // priors
  for (j in 1:A2)
    logitlambdaGP[j] ~ normal(0, 1000);

  // likelihood: main survey
  for (j in 1:N)
    Y[j] ~ bernoulli(1 - exp(-cum_foi[age[j]]));

  // likelihood: Tesh survey
  for (j in 1:NTESH)
    YTESH[j] ~ bernoulli(1 - exp(-meancumfoiG[GTESH[j]]));
}

generated quantities {
  matrix[A, A] cum_foi_by_year;

  // initialize first row and column
  for (i in 1:A) {
    cum_foi_by_year[1, i] = lambda[A - i + 1];
    cum_foi_by_year[i, 1] = lambda[A];
  }

  // fill in the rest
  for (j in 2:A) {
    for (i in 2:A)
      cum_foi_by_year[j, i] = cum_foi_by_year[j-1, i-1] + lambda[A - i + 1];
  }
}
