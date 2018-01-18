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

run_masscyto_model <- function(filename) {
    lvls <- levels(jags_data$c_typ)
    sims <- list()
    for (n in 1:length(lvls)) {
        data_by_cell <- subset(masscyto_model_data, c_typ == lvls[n])
        data_by_cell$conc <- as.numeric(as.character(data_by_cell$conc))
        jags_data <- as.list(data_by_cell)
        n_samps=dim(data_by_cell)[1]
        n_cont_ags <- length(levels(data_by_cell$cont_ag))
        n_exps <- length(levels(data_by_cell$exp))
        jags_data_n <-llist(n_samps, n_cont_ags, n_exps)
        jags_data <- c(jags_data, jags_data_n)
        masscyto_model <- jags.model(file = filename, data = jags_data, n.chains = 3, n.adapt = 1e3)
        update(masscyto_model, 1e3)
        var_names <- c("b0", "b1", "b2")
        masscyto_sim <- coda.samples(model = masscyto_model, n.iter = 9e3, thin=3, 
                                     variable.names = var_names)
        masscyto_csim <- as.mcmc(do.call(rbind, masscyto_sim))
        DIC <- dic.samples(masscyto_model, 3e3, 3)
        Gelman <- gelman.diag(masscyto_sim)
        Raftery <- raftery.diag(masscyto_sim)
        Autocorr <- autocorr.diag(masscyto_sim)
        sims[[n]] <- list(sim=masscyto_sim, csim=masscyto_csim, DIC=DIC, 
                            Gelman=Gelman, Autocorr=Autocorr)
    }
    return(sims)
}
model_filenames <- c("masscyto_model_celltype_dnorm_2.txt", "masscyto_model_celltype_dgamma_2.txt",
                     "masscyto_model_celltype_dlnorm_2.txt")
for (m in 1:length(model_filenames)) {
    models[[m]] <- run_masscyto_model(model_filenames[m]) 
    a <- 1
}

save(list=c("models"), file="masscyto_bycell_modelling_results.Rdata")