require(HDInterval)
c_typs <- levels(masscyto_model_data$c_typ)
exps <- levels(masscyto_model_data$exp)
cont_ags <- levels(masscyto_model_data$cont_ag)
concs <- levels(masscyto_model_data$conc)

concs = c(0.1, 0.3, 1)
predicted <- data.frame(
    contrast_agent = character(),
    cell_type = character(),
    conc = double(),
    pred = double(),
    exp1 = double(),
    exp2 = double(),
    exp3 = double()
)
# gamma
for (conc in concs) {
    for (m in 1:6) {
        model <- models[[2]][[m]]$csim
            col_1_name <- paste0("b0")
            for (b1 in 2:4) {
                col_2_name <- paste0("b1[", b1, "]")
                originals <- mass_cyto_tall[which(
                    mass_cyto_tall$ContrastAgent == cont_ags[b1] &
                        mass_cyto_tall$Concentration == conc &
                        mass_cyto_tall$headers_cell == c_typs[m] &
                        mass_cyto_tall$headers_meas == "Median"
                ),]
                samples <- exp(model[,col_1_name] + model[,col_2_name]*conc)
                prediction <- data.frame(
                    contrast_agent=cont_ags[b1], 
                    cell_type=c_typs[m],
                    conc=conc,
                    pred=median(samples),
                    exp1=originals$value[1],
                    exp2=originals$value[2],
                    exp3=originals$value[3]
                )
                predicted <- rbind(predicted, prediction)
            }
    }
}

pred_tall <- gather(predicted, "meas", "val", 4:7, factor_key=T)

pred_plot <- ggplot(pred_tall) +
    geom_line(aes(x=conc, y=val, col=contrast_agent, 
                  linetype=meas), size=1) +
    facet_grid(.~ cell_type) +
    labs(title="Medians, predicted and by experiment", 
            x="Concentration", y="Signal")

print(pred_plot)