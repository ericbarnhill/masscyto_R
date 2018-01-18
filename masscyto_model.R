# Model 1 - Gamma, random slope for conc
# Model 2 - Normal, random slope for conc
# Model 3 - Gamma, random intercepts only
# Model 4 - Normal, random intercepts only
# Model 5 - Gamma, two random slopes
# Model 6 - Normal, two random slopes
# Model 7 - Gamma, only cell type and cell type * concentration
# Model 8 - Gamma, cell type, cell type * concentration,contrast agent, cont ag * concentration


require(rjags)
require(Hmisc)
jags_data <- as.list(masscyto_model_data)
n_samps=dim(masscyto_model_data)[1]
n_c_typs <- length(levels(masscyto_model_data$c_typ))
n_cont_ags <- length(levels(masscyto_model_data$cont_ag))
n_exps <- length(levels(masscyto_model_data$exp))
jags_data_n <-llist(n_samps, n_c_typs, n_cont_ags, n_exps)
jags_data <- c(jags_data, jags_data_n)
run_masscyto_model <- function(filename) {
  masscyto_model <- jags.model(file = filename, data = jags_data, n.chains = 3, n.adapt = 1e5)
  update(masscyto_model, 1e6)
  var_names <- c("b0", "b1", "b2", "b3")
  masscyto_sim <- coda.samples(model = masscyto_model, n.iter = 1e7,thin=1000, 
                                 variable.names = var_names)
  masscyto_csim <- as.mcmc(do.call(rbind, masscyto_sim))
  DIC <- dic.samples(masscyto_model, 3e3, 3)
  Gelman <- gelman.diag(masscyto_sim)
  Raftery <- raftery.diag(masscyto_sim)
  Autocorr <- autocorr.diag(masscyto_sim)
  model_results <- list(sim=masscyto_sim, DIC=DIC, 
                        Gelman=Gelman, Autocorr=Autocorr)
}
mods <- list()
for (n in 10) {
  model_filename <- paste0("masscyto_model_",n,".txt")
  mods[[n]] <- run_masscyto_model(model_filename) 
}
save(list=c("mods"), file="masscyto_modelling_results.Rdata")