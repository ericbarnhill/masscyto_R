library(EasyABC)
library(gdata)
library(tidyr)
library(rlist)
library(nleqslv)
prep_group_pars <- function(group_stats) {
    mn <- group_stats$value[which(group_stats$MeasurementType == "Mean")]
    mdn <- group_stats$value[which(group_stats$MeasurementType == "Median")]
    pct_05 <- group_stats$value[which(group_stats$MeasurementType == "pct_05")]
    pct_95 <- group_stats$value[which(group_stats$MeasurementType == "pct_95")]
    sum_stats <- c(mn, mdn, pct_05, pct_95)
    return(sum_stats)
}


est_gamma_mle <- function(group_sum_stats) {
    ofn <- function(x,q) {
        sum(abs(c(
            q[c(3,1,4)]-qgamma(c(0.05,0.5,0.95),x[1],x[2]),
            q[2] - median(rgamma(10000, x[1], x[2]))
        )^2))
    }
    dist <- optim(c(1,1),fn=ofn,q=group_sum_stats)
}

gamma_model <- function(par){ 
    samples <- rgamma(10000, shape=par[1], rate=par[2])
    mn <- mean(samples)
    mdn <- median(samples)
    qtls <- quantile(samples, c(0.05, 0.95))
    sum_stats <- c(mn, mdn, qtls[1], qtls[2])
    return(sum_stats)
}

est_gamma_abc <- function(group_sum_stats) {
    ABC_mcmc(method="Marjoram", model=gamma_model, 
             prior=list(c("unif",1e-3,10),c("unif",1e-3,5)), 
             summary_stat_target=group_sum_stats, n_rec = N_REC, 
             verbose=TRUE)
}


generate_distributions <- function(mass_cyto_tall) {
    cell_types <- levels(mass_cyto_tall$CellType)
    n_celltypes <- length(cell_types)
    
    experiments <- levels(mass_cyto_tall$Experiment)
    contrast_agents <- levels(mass_cyto_tall$ContrastAgent)
    concentrations <- levels(mass_cyto_tall$Concentration)
    
    dist_models <- list()
    for (i in 1:length(experiments)) {
        for (j in 1:length(contrast_agents)) {
            for (k in 1:length(concentrations)) {
            experimental_group <- mass_cyto_tall[which(
                mass_cyto_tall$Experiment == experiments[i] &
                mass_cyto_tall$ContrastAgent == contrast_agents[j] &
                mass_cyto_tall$Concentration == concentrations[k]
            ),]
                for (n in 1:n_celltypes) {
                    print(paste("Loop counter: ", i, j, k, n))
                    group_stats <- experimental_group[which(experimental_group$CellType == cell_types[n]),]
                    group_sum_stats <- prep_group_pars(group_stats)
                    model_name <- paste(experiments[i], contrast_agents[j], concentrations[k], cell_types[n])
                    model_key <- paste0(i,j,k,n)
                    if (length(group_sum_stats) == 4) {
                        print(paste("Actual:", group_sum_stats[1], group_sum_stats[2]))
                        N_REC <- 10000
                        m1 <- est_gamma_mle(group_sum_stats)
                        samps1 <- rgamma(10000, m1$par[1], m1$par[2])
                        mn1 <- mean(samps1)
                        mdn1 <- median(samps1)
                        print(paste("Predicted by OPTIM:", format(mn1, digits = 2), format(mdn1, digits = 2)))
                        mod1 <- c(m1, mn1, mdn1)
                        start_time <- Sys.time()
                        m2 <- est_gamma_abc(group_sum_stats)
                        print(paste('Time for ABC calculation:', format(Sys.time()-start_time, digits=2)))
                        params <- m2$param[round(N_REC/2):N_REC,] # leave first half as burn in
                        samps2 <- rgamma(10000, mean(params[,1]), mean(params[,2]))
                        mn2 <- mean(samps2)
                        mdn2<- median(samps2)
                        print(paste("Predicted by ABC:", format(mn2, digits = 2), format(mdn2, digits = 2)))
                        mod2 <- c(m2, mn2, mdn2)
                        mod <- list(mod_mle=mod1, mod_abc=mod2, name=model_name, key=model_key)
                        dist_models <- list.append(dist_models, mod)
                        print("---")
                    }
                }
            }
        }
    }
    return(dist_models)
}


