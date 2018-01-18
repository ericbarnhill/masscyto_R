library(ggplot2)
library(rlist)
# whoops, go get original summary stats
orig_sum_stats <- list()

for (i in 1:length(experiments)) {
    for (j in 1:length(contrast_agents)) {
        for (k in 1:length(concentrations)) {
            experimental_group <- mass_cyto_tall[which(
                mass_cyto_tall$Experiment == experiments[i] &
                mass_cyto_tall$ContrastAgent == contrast_agents[j] &
                mass_cyto_tall$Concentration == concentrations[k]
            ),]
            for (n in 1:n_celltypes) {
                group_stats <- experimental_group[which(experimental_group$headers_cell == cell_types[n]),]
                group_sum_stats <- prep_group_pars(group_stats)
                if (length(group_sum_stats) == 4) {
                    orig_sum_stats <- list.append(orig_sum_stats, group_sum_stats)
                }
            }
        }
    }
}



N <- length(models)
means <- list(or=vector(length=N), m1=vector(length=N), m2=vector(length=N), m3=vector(length=N))
medians <- list(or=vector(length=N), m1=vector(length=N), m2=vector(length=N), m3=vector(length=N))

for (n in 1:N) {
    means$or[n] <- orig_sum_stats[[n]][1]
    means$m1[n] <- models[[n]]$mod1[[9]]
    means$m2[n] <- models[[n]]$mod2[[6]]
    means$m3[n] <- models[[n]]$mod3[[9]]
    
    medians$or[n] <- orig_sum_stats[[n]][2]
    medians$m1[n] <- models[[n]]$mod1[[10]]
    medians$m2[n] <- models[[n]]$mod2[[7]]
    medians$m3[n] <- models[[n]]$mod3[[10]]
    
}

means_df <- data.frame(means)
means_df <- cbind(rep("Mean", dim(means_df)[1]), means_df)
colnames(means_df)[1] <- "meas"
medians_df <- data.frame(medians)
medians_df <- cbind(rep("Median", dim(medians_df)[1]), medians_df)
colnames(medians_df)[1] <- "meas"

means_tall <- gather(means_df, key="model", value="value", -meas)
medians_tall <- gather(medians_df, key="model", value="value",-meas)
all_tall <- rbind(means_tall, medians_tall)
x_vals <- rep(seq(1, N, 1),4*2)

get_resid <- function(x,y) {
    resid <- sqrt(mean(sum((x[,2] - x[,y])^2))) 
}
mean_resids <- vector(length=3)
median_resids <- vector(length=3)

for (n in 1:3) {
    mean_resids[n] <- get_resid(means_df, n+2)
    median_resids[n] <- get_resid(medians_df, n+2)
}

print(paste("Mean RMSE: ", mean_resids))
print(paste("Median RMSE: ", median_resids))


ggplot(all_tall) + 
    geom_line(aes(x=x_vals, y=value, group=model, color=model)) +
    facet_grid(. ~ meas)