load_rtdc_3 <- function(XL_PATH) {
    col_names <- c(
        "cell_type",
        "cont_ag",
        "pooled_controls",
        "subj",
        "prot",
        "mean_area",
        "sd_area", 
        "05_area",
        "95_area", 
        "mean_def",
        "sd_def",
        "05_def",
        "95_def",
        "mean_e",
        "sd_e",
        "05_e",
        "95_e"
    )
    n_sheets = gdata::sheetCount(XL_PATH)
    sheet_names <- gdata::sheetNames(XL_PATH)
    rtdc = list(length=n_sheets-1)
    # skip first sheet -- just contains var names
    for (s in 1:n_sheets) {
        sheet = read.xls(XL_PATH, header=T, sheet = s, stringsAsFactors = F)
        sheetname_split <- unlist(strsplit(sheet_names[s], "_"))
        cell_type <- rep(sheetname_split[2], nrow(sheet))
        cont_ag <- rep(sheetname_split[1], nrow(sheet))
        pooled_controls <- cont_ag
        #fix bad spelling
        cell_type <- unname(sapply(cell_type, function(x) {ifelse(x == "Neuthrophils", "Neutrophils", x)}))
        subj_prot <- t(sapply(sheet[,1], function(x) {
            split_list <- unlist(strsplit(x, "_"))
            subj <- split_list[2]
            prot <- split_list[3]
            subj_prot <- unname(cbind(subj,prot))
        }))
        subj_prot <- subj_prot %>%
            unname %>%
            as.data.frame %>%
            set_colnames(c("subj", "prot"))
        subj <- subj_prot$subj
        prot <- clean_prots(subj_prot$prot)
        pooled_controls[which(
            prot == "Control"
        )] <- "Control"
        sheet <- cbind(cell_type, cont_ag, as.character(pooled_controls), subj, prot, sheet[,2:ncol(sheet)])
        colnames(sheet) <- col_names
        rtdc[[s]] <- sheet
    }
    rtdc <-do.call(rbind, rtdc)
    require(tidyr)
    rtdc_tall <- gather(rtdc, key="measurement", value="value", -(1:5) )
    rtdc_tall$measurement <- factor(rtdc_tall$measurement)
    #bad_entries <- which(rtdc_tall$measurement == "mean_e" & rtdc_tall$value > 2)
    #rtdc_tall$value[bad_entries] <- 
    #    rtdc_tall$value[bad_entries] %>% `/` (100)
    return(rtdc_tall)
}

clean_prots <- function(prots_raw) {
    control_entries <- c("Control", "control", "cntl", "Contrl", "ctl", "Ctrol", "cotrl", "Ctrl", "Ctr", "Cntrol", "cotnrol", 
                         "control-", "Control2")
    #gdcl3_entries <- c("GdCl3")
    #gadovist_entries <- c("Gadovist", "Gadov")
    #dotarem_entries <- c("Dotarem", "Dota")
    #magnevist_entries <- c("Magnevist", "Magnevist_I", "Magnevist_II", "Magnevist_III")
    prots <- sapply(as.character(prots_raw), function(x) {
        x_split_1 <- unlist(strsplit(x, " "))
        x_split_2 <- unlist(strsplit(x_split_1[1], "-"))
        x <- x_split_2[1]
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
    prots_2 <- prots %>%
        as.factor %>%
        unname
    return(prots_2)
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

rtdc_boxplot_pooled_controls_treatments <- function(rtdc_data, stat, title, FUN, log) {
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
        facet_grid(. ~ cell_type) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title)
}

rtdc_boxplot_pooled_controls <- function(rtdc_data, stat, title, FUN, log) {
    if (missing(log)) {
        log = F;
    }
    if (log) {
        rtdc_data$value = log(rtdc_data$value)
    }
    ggplot(subset(rtdc_data, measurement == stat)) +
        geom_jitter(aes(x=pooled_controls, y=value), colour="gray60", width=0.1) + 
        stat_summary(aes(x=pooled_controls, y=value), fun.y=FUN, geom="point", shape=18,
                     size=3, color="red") +
        facet_grid(. ~ cell_type) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title)
}