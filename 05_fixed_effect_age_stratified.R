source("01_load_data.R")
source("02_model_setup.R")

# Fit model on strata with nonzero number of deaths for people 9 or under
data[,  keep_9:=sum(deaths_under_9)>0, by=stratum]
model_9 <- gnm(deaths_under_9 ~ cbt, eliminate = stratum, 
               family = quasipoisson(), data = data, subset=keep_9)

# Fit model on strata with nonzero number of deaths for people between 10 and 19
data[,  keep_10_19:=sum(deaths_10_19)>0, by=stratum]
model_10_19 <- gnm(deaths_10_19 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_10_19)

# Fit model on strata with nonzero number of deaths for people between 20 and 34
data[,  keep_20_34:=sum(deaths_20_34)>0, by=stratum]
model_20_34 <- gnm(deaths_20_34 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_20_34)

# Fit model on strata with nonzero number of deaths for people between 35 and 64
data[,  keep_35_64:=sum(deaths_35_64)>0, by=stratum]
model_35_64 <- gnm(deaths_35_64 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_35_64)

# Fit model on strata with nonzero number of deaths for people 65 or over
data[,  keep_65:=sum(deaths_65_plus)>0, by=stratum]
model_65 <- gnm(deaths_65_plus ~ cbt, eliminate = stratum, 
                family = quasipoisson(), data = data, subset=keep_65)

# Make predictions with all models
pred_9 <- crosspred(cbt, model_9, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_10_19 <- crosspred(cbt, model_10_19, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_20_34 <- crosspred(cbt, model_20_34, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_35_64 <- crosspred(cbt, model_35_64, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_65 <- crosspred(cbt, model_65, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Visually compare results of all models
col <- c("black",'firebrick', "darkorange2",'cadetblue','royalblue')
par(mfrow = c(1, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_9, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.005)))

lines(pred_10_19, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.005)))

lines(pred_20_34, "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.005)))

lines(pred_35_64, "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.005)))

lines(pred_65, "overall", lag=0, cumul=TRUE, ci="bars",col=col[5], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[5], 0.005)))

legend("topleft", c("9 and under","10-19", "20-34", "35-64","65 and above"), lty=c(1,2,3,4,5,6), lwd=1, col=col,
       bty="n", inset=0.05, y.intersp=1, cex=0.7)

hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.15))