data {
  int<lower = 1> N;                  // number of data points
  vector[N] logvalue;                   // dep variable
  int<lower = 1> J;                  // number of cell types
  vector[N] cconc;                   // concentration    
  int<lower = 1, upper = J> typ[N];  // cell id
  vector<lower = 0>[N] SD;           // sd 
}

parameters {
  vector[N] true_value;        // true unknown value 
  vector[2] beta;              //fixed intercept and slope
  real<lower=0> sigma_e;       //error sd
  vector<lower=0>[2] sigma_u;       //cell type sds
  cholesky_factor_corr[2] L_u;
  matrix[2,J] z_u;                     
}

transformed parameters {
  vector[N] mu;
  matrix[2,J] u;
  matrix[J,2] u_t; //transformed u: u'
  u = diag_pre_multiply(sigma_u,L_u)*z_u;  
  u_t = u';
  mu = beta[1] + u_t[typ,1] + cconc .* (beta[2] + u_t[typ,2]);
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
  vector[J] diff_by_celltype;
  matrix[2,2] rho_u;
  rho_u = L_u * L_u';
  for(i in 1:J){
   diff_by_celltype[i] =  (exp(beta[1] + u_t[i,1] + beta[2]+u_t[i,2])) - exp(beta[1] + u_t[i,1]);
  } 
  Conc = (exp(beta[1] + beta[2]) - exp(beta[1]));
}

