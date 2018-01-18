N_SAMPS <- 100
build_normalized_dists <- function(models) {
    masscyto_model_data <- data.frame(
        exp = character(),
        cont_ag= character(),
        conc = character(), 
        c_typ = character(),
        value = double(),
        stringsAsFactors = FALSE
    )
    N <- length(models)
    for (n in 1:N) {
        mod <- models[[n]]
        mod_key <- mod$key
        exp <- experiments[as.numeric(substr(mod_key,1,1))]
        cont_ag <- contrast_agents[as.numeric(substr(mod_key,2,2))]
        conc <- concentrations[as.numeric(substr(mod_key,3,3))]
        c_typ <- cell_types[as.numeric(substr(mod_key,4,4))]
        mod_params <- mod$mod_abc$param[round(N_SAMPS/2):N_SAMPS,]
        shape <- mean(mod_params[,1])
        rate <- mean(mod_params[,2])
        samples <- rgamma(N_SAMPS, shape, rate)
        dist <- data.frame(rep(exp, N_SAMPS), rep(cont_ag, N_SAMPS), rep(conc, N_SAMPS), 
                      rep(c_typ, N_SAMPS), samples)
        colnames(dist) <- c("exp", "cont_ag", "conc", "c_typ", "value")
        masscyto_model_data <- rbind(masscyto_model_data, dist)
    }
    return(masscyto_model_data)
}

masscyto_model_data <- build_normalized_dists(models)
