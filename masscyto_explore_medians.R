require(HDInterval)


explore_medians <- function(mass_cyto_tall, masscyto_model_data, post_models) {
    # if samples data is missing, only print original medians
    if (missing(masscyto_model_data)) {
        masscyto_model_data <- NULL
    }
    c_typs <- levels(mass_cyto_tall$CellType)
    exps <- levels(mass_cyto_tall$Experiment)
    cont_ags <- levels(mass_cyto_tall$ContrastAgent)
    concs <- levels(mass_cyto_tall$Concentration)
    
    concs = c(0.1, 0.3, 1)
    medians <- data.frame(
        contrast_agent = character(),
        cell_type = character(),
        conc = double(),
        pred = double(),
        exp1 = double(),
        exp2 = double(),
        exp3 = double()
    )
    # gamma
    for (m in 1:6) {
        if (!is.null(masscyto_model_data)) {
            # model 2 is the gamma and has lowest DIC
            samples <- post_models[[2]][[m]]$csim
        }
        for (conc in concs) {
               col_1_name <- paste0("b0")
                for (b1 in 2:4) {
                    col_2_name <- paste0("b1[", b1, "]")
                    originals <- mass_cyto_tall[which(
                        mass_cyto_tall$ContrastAgent == cont_ags[b1] &
                            mass_cyto_tall$Concentration == conc &
                            mass_cyto_tall$CellType == c_typs[m] &
                            mass_cyto_tall$MeasurementType == "Median"
                    ),]
                    # inverse of log link
                    if (!is.null(masscyto_model_data)) {
                        samples <- exp(samples[,col_1_name] + samples[,col_2_name]*conc)
                    } 
                    if (!is.null(masscyto_model_data)) {
                        pred=median(samples)
                    } else {
                        pred = 0
                    }
                    prediction <- data.frame(
                        contrast_agent=cont_ags[b1], 
                        cell_type=c_typs[m],
                        conc=conc,
                        pred=pred,
                        exp1=originals$value[1],
                        exp2=originals$value[2],
                        exp3=originals$value[3]
                    )
                    medians <- rbind(medians, prediction)
                }
        }
    }
    
    pred_tall <- gather(medians, "meas", "val", 4:7, factor_key=T)
    # if concentration is treated as factor, need to use points
    #pred_tall$conc <- as.factor(pred_tall$conc)
    
    if (is.null(masscyto_model_data)) {
        pred_tall <- subset(pred_tall, meas != "pred")
        title = "Medians As Reported By Cytometer"
    } else {
        title = "Medians, Reported And medians By Model"
    }
    pred_plot <- ggplot(pred_tall) +
        geom_line(aes(x=conc, y=val, col=contrast_agent, 
                      linetype=meas), size=1) +
        facet_wrap(~ cell_type, scales = "free") +
        labs(title=title, 
                x="Concentration", y="Signal")
    
    print(pred_plot)
    
    return(medians)
}