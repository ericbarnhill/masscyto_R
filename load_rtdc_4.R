load_data <- function(path) {
    col_names <- c(
        "Subject",
        "Date",
        "Cell_Type",
        "Contrast_Agent",
        "Protocol",
        "Events",
        "E_Mean",
        "E_SD", 
        "E_q05",
        "E_q95"
    )
    n_sheets = gdata::sheetCount(path)
    sheet_names <- gdata::sheetNames(path)
    rtdc = list(length=n_sheets-1)
    for (s in 1:n_sheets) {
        sheet = read.xls(path, header=T, sheet = s, stringsAsFactors = F)
        
        # get Cell_Type and Contrast_Agent
        sheetname_split <- unlist(strsplit(sheet_names[s], "_"))
        cell_type <- rep(sheetname_split[2], nrow(sheet))
        cont_ag <- rep(sheetname_split[1], nrow(sheet))
        
        # fix bad spelling
        cell_type <- unname(sapply(cell_type, function(x) {ifelse(x == "Neuthrophils", "Neutrophils", x)}))
        
        # get Date, Subject, Protocol (control or treatment)
        date_subj_prot <- t(sapply(sheet[,1], function(x) {
            split_list <- unlist(strsplit(x, "_"))
            dates <- split_list[1]
            subj <- split_list[2]
            prot <- split_list[3]
            date_subj_prot <- unname(cbind(dates,subj,prot))
        }))
        
        # split subject from protocol
        date_subj_prot <- date_subj_prot %>%
            unname %>%
            as.data.frame %>%
            set_colnames(c("date", "subj", "prot"))
        date <- date_subj_prot$date
        subj <- as.character(date_subj_prot$subj)
        prot <- as.character(clean_prots(date_subj_prot$prot))
        events <- as.numeric(sheet$Events)
        
        # fix two representations for Angela
        subj <- unname(sapply(subj, function(x) {ifelse(x == "AAS", "AngelaS", x)}))
        
        # split Date from colon
        date <- sapply(date, function(x) {
            split_date <- unlist(strsplit(as.character(x), ' '))
            if (length(split_date) > 1) {
                date = split_date[2]
            }
        })
        
        #combine labels with Young's modulus data
        sheet_clean_labels <- data.frame(subj, date, cell_type, cont_ag, prot, events)
        youngs_cols <- unname(sapply(colnames(sheet), function(x) {grepl("Young", x)}))
        youngs_data <- sheet[,youngs_cols]
        rtdc_data <- cbind(sheet_clean_labels, youngs_data)
        colnames(rtdc_data) <- col_names
        
        #calculate conservativce SE
        rtdc_data <- mutate(rtdc_data, E_SD_hi_log = log(E_q95) - log(E_Mean)) %>%
            mutate(., E_SD_low_log = log(E_Mean) - log(E_q05)) %>%
            mutate(., E_SE_conserv_log = E_SD_hi_log / sqrt(Events)) %>%
            mutate(., E_Mean_log = log(E_Mean))
            
        
        rtdc[[s]] = rtdc_data
    }
    rtdc <-do.call(rbind, rtdc)
    cols_to_factor = colnames(rtdc)[1:5]
    rtdc[cols_to_factor] <- lapply(rtdc[cols_to_factor], function(x) {factor(unlist(x))})
    rownames(rtdc) <- c()
    return(rtdc)    
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

