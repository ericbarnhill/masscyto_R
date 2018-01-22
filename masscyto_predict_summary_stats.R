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

exp_mean <- vector(length=42)
exp_median <- vector(length=42)
exp_05 <- vector(length=42)
exp_95 <- vector(length=42)

mom_mean <- vector(length=42)
mom_median <- vector(length=42)
mom_05 <- vector(length=42)
mom_95 <- vector(length=42)

mle_mean <- vector(length=42)
mle_median <- vector(length=42)
mle_05 <- vector(length=42)
mle_95 <- vector(length=42)

abc_mean <- vector(length=42)
abc_median <- vector(length=42)
abc_05 <- vector(length=42)
abc_95 <-vector(length=42)
abc_uncertainty.shape <- vector(length=42)
abc_uncertainty.rate <- vector(length=42)


collect_dist_results <- function(mass_cyto_tall, dist_models) {
    N_REC = 10000
    
    cell_types <- levels(mass_cyto_tall$CellType)
    n_celltypes <- length(cell_types)
    
    experiments <- levels(mass_cyto_tall$Experiment)
    contrast_agents <- levels(mass_cyto_tall$ContrastAgent)
    concentrations <- levels(mass_cyto_tall$Concentration)
    
    index = 1;
    for (i in 1:length(experiments)) {
    #for (i in 1:1) { # only analyze first experiment
        for (j in 1:length(contrast_agents)) {
            for (k in 1:length(concentrations)) {
            experimental_group <- mass_cyto_tall[which(
                mass_cyto_tall$Experiment == experiments[i] &
                mass_cyto_tall$ContrastAgent == contrast_agents[j] &
                mass_cyto_tall$Concentration == concentrations[k]
            ),]
                for (n in 1:n_celltypes) {
                    group_stats <- experimental_group[which(experimental_group$CellType == cell_types[n]),]
                    group_sum_stats <- prep_group_pars(group_stats)
                    model_name <- paste(experiments[i], contrast_agents[j], concentrations[k], cell_types[n])
                    model_key <- paste0(i,j,k,n)
                    if (length(group_sum_stats) == 4 & index <= 42) { # kludge till data is re-run
                    #if (length(group_sum_stats) == 4) {
                            #exp
                        exp_mean[index] <- group_sum_stats[1]
                        exp_median[index] <- group_sum_stats[2]
                        exp_05[index] <- group_sum_stats[3]
                        exp_95[index] <- group_sum_stats[4]
                        #mom
                        mom <- dist_models[[index]]$mod1
                        samps <- rgamma(10000, mom$shape, mom$rate)
                        mom_mean[index] <- mean(samps)
                        mom_median[index] <- median(samps)
                        mom_05[index] <- quantile(samps, 0.05)
                        mom_95[index] <- quantile(samps, 0.95)
                        #mle
                        mle <- dist_models[[index]]$mod2
                        samps <- rgamma(1000, mle$par[1], mle$par[2])
                        mle_mean[index] <- mean(samps)
                        mle_median[index] <- median(samps)
                        mle_05[index] <- quantile(samps, 0.05)
                        mle_95[index] <- quantile(samps, 0.95)
                        #
                        abc <- dist_models[[index]]$mod3
                        params <- abc$param[round(N_REC/2):N_REC,] # leave first half as burn in
                        shape <- params[,1]
                        rate <- params[,2]
                        samps <- rgamma(10000, mean(shape), mean(rate))
                        abc_mean[index] <- mean(samps)
                        abc_median[index] <- median(samps)
                        abc_05[index] <-quantile(samps, 0.05)
                        abc_95[index] <-quantile(samps, 0.95)
                        abc_shape_hdi <- hdi(shape)
                        abc_rate_hdi <- hdi(rate)
                        abc_uncertainty.shape[index]  <- 
                            (abc_shape_hdi[2] - abc_shape_hdi[1]) / mean(shape)
                        abc_uncertainty.rate[index]  <- 
                            (abc_rate_hdi[2] - abc_rate_hdi[1]) / mean(rate)
                        index = index + 1
                    }
                }
            }
        }
    }
    exp_cond <- seq(1,length(mom_mean),1)
    dist_results_wide <- data.frame(
        exp_mean, exp_median, exp_05, exp_95,
        mom_mean, mom_median, mom_05, mom_95,
        mle_mean, mle_median, mle_05, mle_95, 
        abc_mean, abc_median, abc_05, abc_95,
        abc_uncertainty.shape, abc_uncertainty.rate,
        exp_cond
    )
    
    error_colnames <- c("mom_mean.err", "mom_median.err", "mom_05.err", "mom_95.err",
                        "mle_mean.err", "mle_median.err", "mle_05.err", "mle_95.err",
                        "abc_mean.err", "abc_median.err", "abc_05.err", "abc_95.err") 
    
    error_results <- matrix(nrow=nrow(dist_results_wide), ncol=length(error_colnames))
    for (i in 1:3) {
        for (j in 1:4) {
            col_num  = (i)*4+j
            #error <- (dist_results_wide[,col_num] - dist_results_wide[,j])^2
            error <- sqrt((dist_results_wide[,col_num] - dist_results_wide[,j])^2) / dist_results_wide[,j]
            error_results[,col_num-4] <- error
        }
    }
    colnames(error_results) <- error_colnames
    
    dist_results <- gather(dist_results_wide, "stat", "value", -exp_cond, factor_key = T)
    dist_results <- separate(data = dist_results, col = "stat",
                             into = c("method", "stat"), sep = "_")
    
    results_list <- list(error_results = error_results, dist_results = dist_results,
                         dist_results_wide = dist_results_wide)
    
}

plot_dist_results <- function(dist_results) {
    require('RColorBrewer')
    cols <- brewer.pal(3, 'Set1')
    ggplot(dist_results) + 
            geom_line(aes(x=exp_cond, y=value, group=method, col=method, linetype = method)) +
        facet_wrap(~ stat, scales="free") +
        scale_linetype_manual(values=c("dashed", "solid", "dashed", "dashed")) +
        labs(title="Comparison of Experimental and Predicted Summary Statistics", 
             x="All Experiments", y="Signal") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

make_error_table <-function(error_results) {
    error_sums <- colMeans(error_results)
    err_melt <- melt(error_sums)
    err_df <- data.frame(rownames(err_melt), err_melt$value)
    err_sep <- separate(err_df, col = 1, into = c("method", "stat"), sep="_")
    err_table <- spread(err_sep, key = method, value=err_melt.value)
}