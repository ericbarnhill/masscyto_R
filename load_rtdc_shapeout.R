load_rtdc <- function(XL_PATH) {
    require(magrittr, dplyr, tidyr)
    col_names <- c(
        "contrast_type",
        "cell_type", 
        "condition",
        "subject",
        "flow rate",
        "% gated",
        "events",
        "median area",
        "mode area",
        "mean area",
        "sd area",
        "median deformation",
        "mode deformation",
        "mean deformation",
        "sd deformation",
        "median E",
        "mode E",
        "mean E",
        "sd E"
    )
    n_sheets = gdata::sheetCount(XL_PATH)
    sheet_names <- gdata::sheetNames(XL_PATH)
    rtdc = list(length=n_sheets-1)
    # skip first sheet -- just contains var names
    num_data_sets_by_sheet <- c(6, 6, 8, 8)
    for (s in 1:n_sheets) {
        sheet_xl = read.xls(XL_PATH, header=F, sheet = s, stringsAsFactors = F)
        contrast_type <- sheet_names[s]
        n_sets = num_data_sets_by_sheet[s]
        sheet = subset(sheet_xl, V4 != '' & !is.na(as.numeric(V4)))
        sheet <- cbind(c(rep("all", n_sets*2), rep("monos", n_sets*2), rep("neutros", n_sets*2)), sheet)
        sheet <- cbind(rep(contrast_type, nrow(sheet)), sheet)
        sheet[,3] <- rep(c(rep("treatment", n_sets), rep("control", n_sets)), 3)
        sheet[,-c(1:4)] <- lapply(sheet[,-c(1:4)], as.numeric)
        colnames(sheet) <- col_names
        print("")
    }
    rtdc <-do.call(rbind, rtdc)
    require(tidyr)
    rtdc_tall <- gather(rtdc, key="measurement", value="value", -(1:4) )
    rtdc_tall$measurement <- factor(rtdc_tall$measurement)
    bad_entries <- which(rtdc_tall$measurement == "mean_e" & rtdc_tall$value > 2)
    rtdc_tall$value[bad_entries] <- 
        rtdc_tall$value[bad_entries] %>% `/` (100)
    return(rtdc_tall)
}

clean_prots <- function(prots) {
    control_entries <- c("cntl", "Contrl", "ctl", "Ctrol", "cotrl", "Ctrl", "Ctr", "Cntrol")
    #gdcl3_entries <- c("GdCl3")
    #gadovist_entries <- c("Gadovist", "Gadov")
    #dotarem_entries <- c("Dotarem", "Dota")
    #magnevist_entries <- c("Magnevist", "Magnevist_I", "Magnevist_II", "Magnevist_III")
    prots <- sapply(as.character(prots), function(x) {
        if (x %in% control_entries) {
            x <- "Control"
        } else {
            x <- "Treatment"
        }
        #} else if (x %in% gadovist_entries) {
        #    x <- "Gadovist"
        #} else if (x %in% dotarem_entries) {
        #    x <- "Dotarem"
        #} else if (x %in% magnevist_entries) {
        #    x <- "Magnevist"
        #} else if (x %in% gdcl3_entries) {
        #    x <- "GdCl3"
        #} else {
        #    print("ERROR: entry not labeled properly as control or contrast agent")
        #}
        return(x)
    })
    prots %>%
        as.factor %>%
        unname
}

rtdc_exploratory_boxplot <- function(rtdc_data, stat, title, FUN, log) {
    if (missing(log)) {
        log = F;
    }
    if (log) {
        rtdc_data$value = log(rtdc_data$value)
    }
    ggplot(subset(rtdc_data, measurement == stat)) +
        geom_jitter(aes(x=prot, y=value), colour="gray60", width=0.1) + 
        stat_summary(aes(x=prot, y=value), fun.y=FUN, geom="point", shape=18,
                     size=3, color="red") +
        facet_grid(. ~ cell_type + cont_ag) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title)
}