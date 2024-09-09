source("01_load_data.R")
source("02_model_setup.R")

# Fit model on strata with nonzero number of car deaths
data[,  keep_veh:=sum(median_vehicle_round)>0, by=stratum]
model_veh <- gnm(median_vehicle_round ~ cbt, eliminate = stratum, 
                 family = quasipoisson(), data = data, subset=keep_veh)

# Fit model on strata with nonzero number of motorcycle deaths
data[,  keep_mot:=sum(median_motorcycle_round)>0, by=stratum]
model_mot <- gnm(median_motorcycle_round ~ cbt, eliminate = stratum, 
                 family = quasipoisson(), data = data, subset=keep_mot)

# Fit model on strata with nonzero number of bike deaths
data[,  keep_bike:=sum(median_bicycle_round)>0, by=stratum]
model_bike <- gnm(median_bicycle_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_bike)

# Fit model on strata with nonzero number of pedestrian deaths
data[,  keep_ped:=sum(median_ped_round)>0, by=stratum]
model_ped <- gnm(median_ped_round ~ cbt, eliminate = stratum, 
                 family = quasipoisson(), data = data, subset=keep_ped)

# Make predictions with all models
pred_veh <- crosspred(cbt, model_veh, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_mot <- crosspred(cbt, model_mot, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_bike <- crosspred(cbt, model_bike, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_ped <- crosspred(cbt, model_ped, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Visually compare results of all models
col <- c('firebrick', "darkorange2",'cadetblue','royalblue')
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_veh, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))

lines(pred_mot, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))

lines(pred_bike, "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))

lines(pred_ped, "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))

legend("topleft", c("car","motorcycle", "bike", "pedestrian"), lty=c(1,2,3,4,5,6), lwd=1, col=col,
       bty="n", inset=0.05, y.intersp=1, cex=0.7)

hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.15))
