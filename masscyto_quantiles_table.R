quantiles_table <- function(mass_cyto_tall) {
    c_typs <- levels(mass_cyto_tall$CellType)
    exps <- levels(mass_cyto_tall$Experiment)
    cont_ags <- levels(mass_cyto_tall$ContrastAgent)
    concs <- levels(mass_cyto_tall$Concentration)
    measurements = levels(mass_cyto_tall$MeasurementType)[2:5]
    quantiles <- matrix(nrow=0,ncol=7)
    concs <- c(0.1, 0.3, 1.0)
    for (m in 1:6) {
        for (c in 1:length(concs)) {
            conc <- concs[c]
            col_1_name <- paste0("b0")
            for (b1 in 2:4) {
                col_2_name <- paste0("b1[", b1, "]")
                condition <- mass_cyto_tall[which(
                    mass_cyto_tall$ContrastAgent == cont_ags[b1] &
                        mass_cyto_tall$Concentration == conc &
                        mass_cyto_tall$CellType == c_typs[m]
                ),]
                meas_vec <- vector(length=6)
                for (n in 1:length(measurements)) {
                    meas_vec[n] = mean(condition[which(
                        condition$MeasurementType == measurements[n]
                    ),"value"])
                }
                meas_vec[5] <- c_typs[m]
                meas_vec[6] <- b1-1
                meas_vec[7] <- conc
                quantiles <- rbind(quantiles, meas_vec)
            }
        }
    }
    quantiles <- as.data.frame(quantiles)
    colnames(quantiles) <- c(measurements, "CellType", "Experiment", "Concentration")
    rownames(quantiles) <- NULL
    for (n in c(1:4)) {
        quantiles[,n] <- as.numeric(as.character(quantiles[,n]))
    }
    quantiles_tall <- gather(quantiles, "value_type", "value", factor_key = T, -(5:7))
    #combine the three experiments
    quantiles_tall <- aggregate(value ~ CellType + Concentration + value_type, quantiles_tall, mean)
    quantiles_plot <- ggplot(quantiles_tall) +
        geom_line(aes(x=Concentration, y=value, group=value_type, col=value_type)) +
        labs(title="Cytometer Quantiles for Experimental Conditions", 
             x="Concentration x Experiment", y="Signal") +
        facet_wrap(~ CellType, scales="free", ncol=3) +
        scale_x_discrete(expand = c(0,0))
    print(quantiles_plot)
    return(quantiles_tall)
}
