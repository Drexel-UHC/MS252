source("01_load_data.R")

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
# argvartmean <- list(fun = "ns", df = 3)
# knotstmean <- equalknots(data[[Temp_measure]], df=3)
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)

# Create crossbasis for DLNM
cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                  argvar = argvartmean, arglag = list(knots = lagknots))

# Specify stratification scheme
# stratify by city-year-month-dow
data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]

# Specify which rows to keep
# keep all strata with nonzero number of road deaths
data[,  keep:=sum(median_road_round)>0, by=stratum]

# Fit model
model_road <- gnm(median_road_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)
# Estimate association via the DLNM basis
pred_road <- crosspred(cbt, model_road,
                       at = seq(minT, maxT, by = .1),
                       cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Plot overal RR and temperature histogram
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road, "overall", lag=0, cumul=TRUE,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road$predvar<=medT
ind2 <- pred_road$predvar>=medT
lines(pred_road$predvar[ind1],pred_road$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road$predvar[ind2],pred_road$allRRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))

#y_axis_labels <- seq(0, 1.0, by = 0.05)
#axis(2, at = y_axis_labels)  # Convert to percentage


# dev.off()
# plot(pred_road, xlab="Temperature", zlab="RR", theta=40, phi=25, lphi=30, border=0, shade=3, 
#      main="Heat effects on road mortality frequency", zlim=c(0.9,1.2))
# fqaic(model_road)

# pred_road$allRRfit[c('1','25','75','90','99')]
# pred_road$allRRhigh[c('1','25','75','90','99')]
# pred_road$allRRlow[c('1','25','75','90','99')]
