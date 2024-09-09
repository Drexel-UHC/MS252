# source("01_load_data.R")

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