data {
  int <lower=0> A; //the number of age classes
  int <lower=0> A2; //the number of foi groups
  int <lower=0> N; //the number of individuals
  int <lower=0> age[N]; // Age
  int <lower=0, upper=1> Y[N]; // Outcome
  int <lower=1, upper=A2> ind_by_age[A]; //
  int <lower=0> NTESH; // Number of individuals in Tesh
  int <lower=0> NGTESH; //Number of age groups in Tesh
  int <lower=0> SizeGTESH[NGTESH]; //Size (in years) of each age groups 
  int <lower=0> GTESH[NTESH]; //Age group of each individual
  int <lower=0, upper=1> YTESH[NTESH]; //Outcome TESH
  int <lower=0> GageTESH[A]; //Age group of each FOI
  int <lower=0> Y1; // year 1
  int <lower=0> YearsFixed; // Number of years with fixed lambda
}

parameters {
  real logitlambdaGP[A2]; //
}

transformed parameters {

  real lambda[A];
  real cum_foi[A];   // cumulative foi by age
  real cum_foi2[A];   // cumulative foi by age for Tesh lot
  real logitlambda[A];
  real logitlambda2[A];
  real meancumfoiG[NGTESH]; // mean FOI by age group
  real sumcumfoiG[NGTESH]; // sum FOI by age group

  //Total FOI by year
  for (j in 1:A) {
    if(logitlambdaGP[ind_by_age[j]]<(-6)){
      logitlambda[j] <- (-6);
    }else{
      logitlambda[j] <- logitlambdaGP[ind_by_age[j]];
  }
}
 

  for (k in 1:A) {
        if ((A-YearsFixed)<k){
      logitlambda2[k] <- logitlambda[(A-YearsFixed)];
    }else{
      logitlambda2[k] <- logitlambda[k];
    }
  }

  for (j in 1:A) {
      lambda[j] <- inv_logit(logitlambda2[j]);
    }

  // lambda[1] <- 0;
  // lambda[2] <- 0;
  // lambda[3] <- 0;
  // lambda[4] <- 0;
  // lambda[5] <- 0;

  cum_foi[1] <- lambda[1];
  for (j in 2:A) {
      cum_foi[j] <- cum_foi[j-1]+lambda[j];
    }

  for (j in 1:A) {
    if(j<Y1){
      cum_foi2[j] <- 0;
    }else{
      cum_foi2[j] <- cum_foi[j]-cum_foi[Y1];
    } 
  }

  for (j in 1:NGTESH){
    sumcumfoiG[j] <- 0;
  }

  for (j in 1:A){
    if(GageTESH[j]!=0){
      sumcumfoiG[GageTESH[j]] <- sumcumfoiG[GageTESH[j]]+cum_foi2[j];
    }
  }

  for (j in 1:NGTESH){
    meancumfoiG[j] <- sumcumfoiG[j]/SizeGTESH[j];
  }
}


model {

  //FOI by group
  for (j in 1:A2) {
    logitlambdaGP[j] ~ normal(0, 1000);
  }


  for (j in 1:N) {  
    Y[j] ~ bernoulli(1-exp(-cum_foi[age[j]]));
  }

  for (j in 1:NTESH) {  
    YTESH[j] ~ bernoulli(1-exp(-meancumfoiG[GTESH[j]]));
  }

  // for (j in 1:NTESH) {
  //   print(GTESH[j]);
  // }

  // print(logitlambda[50]);
  // print(lambda[50]);
  // print(cum_foi[50]);
  // print(meancumfoiG[1]);

}
 



generated quantities {

  matrix[A,A] cum_foi_by_year;

  for (i in 1:A){ 
      cum_foi_by_year[1,i] <- lambda[A-i+1];
      cum_foi_by_year[i,1] <- lambda[A];  //Give everbody in year 1 the same foi
    }

  for (j in 2:A){
    for (i in 2:A){ 
      cum_foi_by_year[j,i] <- cum_foi_by_year[j-1,i-1]+lambda[A-i+1];
    }
   } 


}


