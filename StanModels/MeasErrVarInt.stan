data {
  int<lower = 1> N;                  // number of data points
  vector[N] logvalue;                // dep variable
  int<lower = 1> J;                  // number of types
  vector[N] cconc;                   // concentration    
  vector<lower = 0>[N] SD;           // sd of log value 
  int<lower = 1, upper = J> typ[N];  // type
}

parameters {
  vector[2] beta;              //intercept and slope
  vector[N] true_value;        //true unknown value
  vector[J] u;                 //type varying intercepts
  real<lower=0> sigma_e;       //residual error sd
  real<lower=0> sigma_u;       //type sd
}

model {
  vector[N] mu;
  //priors
  true_value ~ normal(0, 10); 
  beta[1] ~ normal(0, 10);
  beta[2] ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  u ~ normal(0, sigma_u);    //type varying intercepts
  // likelihood
  mu = beta[1] + u[typ] + beta[2] * cconc;
  logvalue ~ normal(true_value, SD); // measurement error
  true_value ~ normal(mu, sigma_e);
}

