require(HDInterval)

contrast_values <- data.frame(
    cell_type_indices = integer(),
    cont_indices = integer(),
    low = double(),
    mean = double(),
    high = double()
)

for (m in 1:6) {
    model_csim <- as.data.frame(models[[2]][[m]]$csim)
    #c21 <- model_csim[["b1[2]"]] - model_csim[["b1[1]"]]
    #c31 <- model_csim[["b1[3]"]] - model_csim[["b1[1]"]]
    #c41 <- model_csim[["b1[4]"]] - model_csim[["b1[1]"]]
    cGadDot <- model_csim[["b1[3]"]] - model_csim[["b1[2]"]]
    cMagDot <- model_csim[["b1[4]"]] - model_csim[["b1[2]"]]
    cMagGad <- model_csim[["b1[4]"]] - model_csim[["b1[3]"]]
    #contrasts <- data.frame(c21, c31, c41, c32, c42, c43)
    contrasts <- data.frame(cGadDot, cMagDot, cMagGad)
    low <- vector(length=ncol(contrasts))
    high <- vector(length=ncol(contrasts))
    mean <- vector(length=ncol(contrasts))
    cont_indices = vector(length=length(contrasts))
    for  (n in 1:ncol(contrasts)) {
        contrast <- contrasts[,n]
        contrast_mean <- mean(contrast)
        contrast_hdi <- hdi(contrast)
        #convert to multiplicate factor
        low[n] <- exp(contrast_hdi[1]) - 1
        high[n] <- exp(contrast_hdi[2]) - 1
        mean[n] <- exp(contrast_mean) - 1
        cont_indices[n] <- n
    }
    cont_vals <- data.frame(cont_indices, low, mean, high)
    cell_type_indices <- rep(m, nrow(cont_vals))
    cont_vals <- cbind(cell_type_indices, cont_vals)
    contrast_values <- rbind(contrast_values, cont_vals)
}

cell_types <- levels(as.factor(headers_cell))
contrast_values$cell_type_indices <- factor(sapply(contrast_values$cell_type_indices, FUN= function(x) cell_types[x]))

plot_obj <- ggplot(contrast_values) + 
    geom_point(aes(x=factor(cont_indices), y=mean, group=cell_type_indices))  + 
    geom_errorbar(aes(x=factor(cont_indices), ymin=low, ymax=high, group=cell_type_indices), alpha = 0.3)  + 
    facet_grid(. ~ cell_type_indices) +
    scale_y_continuous(limits=c(-0.5,1)) +
    scale_x_discrete(labels=factor(colnames(contrasts))) +
    labs(title="Ratios of Contrast Agent Slopes", x="Contrast", y="Ratio")
    

print(plot_obj)