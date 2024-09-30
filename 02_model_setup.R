library("tidyverse")
library("data.table")
library("dlnm"); library("gnm"); library("splines") # Modeling tools

# Load data
path_data <- "C:/Users/dpw48/OneDrive - Drexel University/git/MS252/Data/"

all <- fread(paste0(path_data, "final_temp_cluster_subgroup.csv"))

# Select all temperature clusters for analysis
data <- all

# Select only specific temperature clusters for analysis
# data <- all[all[['cluster_ward_std_6']] %in% c(1)]

# Specify temperature measure
Temp_measure <- "tmp_pw_percentile" # percentile
# Temp_measure <- "L1ADtemp_pw"     # C°
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

# Specify knot number for lags
lagknots <- logknots(3, df = 3)

# Specify knot locations for temperature
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)

# Create crossbasis for DLNM
cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                  argvar = argvartmean, arglag = list(knots = lagknots))

# Specify stratification scheme
# stratify by city-year-month-dow
data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]