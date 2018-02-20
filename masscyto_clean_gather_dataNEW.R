library(gdata)
library(tidyr)
library(rlist)
# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
#data =  rgamma(, mean =5.3, sd = 2.7)

# set some globals
#MASSCYTO_DIR = "/home/ericbarnhill/Documents/code/R/masscyto_R/"
mass_cyto_xls_path <- file.path("experimental_data.xlsx")

clean_measurement_names <- function(mass_cyto_tall) {
    names <- as.character(mass_cyto_tall$MeasurementType)
    is_Median <- sapply(mass_cyto_tall$MeasurementType, function(x) grepl(pattern = "Median", x = x))
    is_Mean <- sapply(mass_cyto_tall$MeasurementType, function(x) grepl(pattern = "Mean", x = x))
    is_95th <- sapply(mass_cyto_tall$MeasurementType, function(x) grepl(pattern = "95th", x = x))
    is_5th <- sapply(mass_cyto_tall$MeasurementType, function(x) grepl(pattern = "5th", x = x))
    is_5th <- is_5th & !is_95th
    names[is_Median] <- "Median"
    names[is_Mean] <- "Mean"
    names[is_95th] <- "pct_95"
    names[is_5th] <- "pct_05"
    mass_cyto_tall$MeasurementType <- factor(names)
    return(mass_cyto_tall)
}

clean_gather_data <- function(mass_cyto) {
    mass_cyto_headers <- read.xls(xls=mass_cyto_xls_path, header=FALSE, skip=1, nrows=2)[,-(5:10)]
    mass_cyto_headers <- apply(X=mass_cyto_headers, MARGIN = 2, FUN = 
                                   function(x) paste(as.character(x[1]), as.character(x[2]), sep="__"))
    colnames(mass_cyto) <- mass_cyto_headers
    colnames(mass_cyto)[1:4] <- c("Sample", "Experiment", "ContrastAgent", "Concentration")
    mass_cyto_tall <- gather(mass_cyto, measurement, value, -Sample, -Experiment, -ContrastAgent, -Concentration)
    mass_cyto_tall$Experiment <- as.character(mass_cyto_tall$Experiment)
    mass_cyto_tall$Experiment[which(mass_cyto_tall$Experiment=="repeat 2")] <- "2"
    mass_cyto_tall$Experiment <- factor(mass_cyto_tall$Experiment)
    mass_cyto_tall$Concentration <- factor(mass_cyto_tall$Concentration)
    headers_split <- unlist(strsplit(mass_cyto_tall$measurement, '__'))
    CellType <- headers_split[seq(1,length(headers_split),2)]
    MeasurementType <- headers_split[seq(2,length(headers_split),2)]
    mass_cyto_tall <- cbind(mass_cyto_tall, CellType, MeasurementType)
    mass_cyto_tall<- subset(mass_cyto_tall, select = -c(measurement, Sample))
    # observed summary statistics
    mass_cyto_tall <- clean_measurement_names(mass_cyto_tall)
    return(mass_cyto_tall)
}

load_clean_data <- function() {
    MASSCYTO_DIR = "/home/ericbarnhill/Documents/code/R/2017-12-15-masscyto/"
    mass_cyto_xls_path <- file.path(MASSCYTO_DIR, "angela.xlsx")
    mass_cyto <- read.xls(xls=mass_cyto_xls_path, header=FALSE, skip=3, nrows=30)[,-(5:10)]
    mass_cyto_tall <- clean_gather_data(mass_cyto)
}