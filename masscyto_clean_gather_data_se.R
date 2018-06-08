library(gdata)
library(tidyr)
library(rlist)
# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
#data =  rgamma(, mean =5.3, sd = 2.7)

# set some globals
#MASSCYTO_DIR = "/home/ericbarnhill/Documents/code/R/masscyto_R/"
#mass_cyto_xls_path <- file.path("experimental_data.xlsx")

set_columns <- function(df, main_columns, stat_types, cell_types, marker_types) {
    #
    n_stat_types <- length(stat_types)
    n_cell_types <- length(cell_types)
    n_marker_types <- length(marker_types)
    stat_types <- rep(stat_types, each = n_cell_types*n_marker_types)
    cell_types <- rep(cell_types, times = n_cell_types, each = n_marker_types)
    marker_types <- rep(marker_types, times = n_stat_types*n_cell_types)
    total_cols <- n_stat_types*n_cell_types*n_marker_types
    meas_names <- unlist(lapply(seq(1, total_cols, 1), 
                         function(x) {
                             paste0(stat_types[x], '__', cell_types[x], '__', marker_types[x])
                         }))
    return(c(main_columns, meas_names))
}

process_events <- function(mass_cyto) {
    cnames <- colnames(mass_cyto)
    events_cols <- unname(sapply(colnames(mass_cyto), function(x) {startsWith(str=x, pattern="Events")}))
    events <- mass_cyto[,events_cols]
    events <- apply(X = events, MARGIN=2, FUN = function(x) {
      x <- x / 100 * mass_cyto["Counts"]
    })
    mass_cyto[,events_cols] <- as.data.frame(events)
    return(mass_cyto)
}


clean_gather_data_se <- function(mass_cyto, main_columns, stat_types, cell_types, marker_types) {
    mass_cyto <- separate(mass_cyto, Sample, c("Subject", "Date"), 
                          sep="_", remove=T, extra="drop", fill="left")
    #
    mass_cyto <- mass_cyto[!is.na(mass_cyto$Subject),]
    mass_cyto_colnames <- set_columns(mass_cyto, main_columns, stat_types, cell_types, marker_types)
    colnames(mass_cyto) <- mass_cyto_colnames
    mass_cyto[,6:ncol(mass_cyto)] <- apply(X = mass_cyto[,6:ncol(mass_cyto)], MARGIN=2, FUN = function(x) {
      as.numeric(as.character(unlist(x)))
    })
    mass_cyto <- process_events(mass_cyto)
    mass_cyto_tall <- gather(mass_cyto, Measurement, Value, -Subject, -Date, -Experiment, -Contrast_Agent, -Concentration, -Counts) %>%
        select(-Counts)
    mass_cyto_tall <- separate(mass_cyto_tall, col = "Measurement", c("Stat_Type", "Cell_Type", "Marker_Type"), sep = "__", remove = T)
    factor_cols <- c("Experiment", "Contrast_Agent", "Stat_Type", "Cell_Type", "Marker_Type")
    mass_cyto_tall[factor_cols] <- lapply(mass_cyto_tall[factor_cols], factor)
    mass_cyto_wide <- spread(mass_cyto_tall, Stat_Type, Value)
    mass_cyto_wide <- mass_cyto_wide %>%
        mutate(., logConcentration = log(Concentration)) %>%
        mutate(., logMean = log(Mean)) %>%
        mutate(., SD_conserv = log(Q95) - log(Mean)) %>%
        mutate(., SE = SD_conserv / sqrt(Events))
    data_list = list(mass_cyto_tall = mass_cyto_tall, mass_cyto_wide = mass_cyto_wide)
    return(data_list)
}

load_clean_data <- function(path) {
    mass_cyto <- read.xls(xls=path)
    main_columns <- c("Subject", "Date", "Experiment", "Contrast_Agent", "Concentration", "Counts")
    stat_types <- c("Events", "Mean", "Median", "CV", "Q05", "Q95")
    cell_types <- c("Neutrophils","T cells","Monocytes","B cells","NK cells_1","NK cells_2")
    marker_types <- c("GdCl3")
    data_list <- clean_gather_data_se(mass_cyto, main_columns, stat_types, cell_types, marker_types)
}