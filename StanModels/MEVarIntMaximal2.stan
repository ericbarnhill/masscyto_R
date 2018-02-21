data {
  int<lower = 1> N;                  // number of data points
  vector[N] logvalue;                   // dep variable
  int<lower = 1> J;                  // number of cell types
  vector[N] cconc;                   // concentration
  vector[N] MvsD;
  vector[N] DvsG;
  vector[N] cconcMvsD;
  vector[N] cconcDvsG;   
  int<lower = 1, upper = J> typ[N];  // cell id
  vector<lower = 0>[N] SD;           // sd 
}

parameters {
  vector[N] true_value;        // true unknown value 
  vector[6] beta;              //fixed intercept and slope
  real<lower=0> sigma_e;       //error sd
  vector<lower=0>[6] sigma_u;       //cell type sds
  cholesky_factor_corr[6] L_u;
  matrix[6,J] z_u;                     
}

transformed parameters {
  vector[N] mu;
  matrix[6,J] u;
  matrix[J,6] u_t; //transformed u: u'
  u = diag_pre_multiply(sigma_u,L_u)*z_u;  
  u_t = u';
  mu = beta[1] + u_t[typ,1] + cconc .* (beta[2] + u_t[typ,2]) + 
  MvsD .* (beta[3] + u_t[typ,3]) + DvsG .* (beta[4] + u_t[typ,4]) +
  cconcMvsD .* (beta[5] + u_t[typ,5]) + cconcDvsG .* (beta[6] + u_t[typ,6]);
}

model {
  //priors
  true_value ~ normal(0, 10);         
  logvalue ~ normal(true_value, SD); // measurement error
  beta[1] ~ normal(0, 10);
  beta[2] ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  to_vector(z_u) ~ normal(0,1);
  L_u ~ lkj_corr_cholesky(2);
  // likelihood
  true_value ~ normal(mu, sigma_e);
}
generated quantities {
  real Conc;
  matrix[6,6] rho_u;
  rho_u = L_u * L_u'; 
  Conc = (exp(beta[1] + beta[2]) - exp(beta[1]));
  // can generate others if needed
}

