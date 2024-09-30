source("01_load_data.R")
source("02_model_setup.R")

# Fit model on strata with nonzero number of male deaths
data[,  keep_male:=sum(male_deaths)>0, by=stratum]
model_male <- gnm(male_deaths ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_male)

# Fit model on strata with nonzero number of female deaths
data[,  keep_female:=sum(female_deaths)>0, by=stratum]
model_female <- gnm(female_deaths ~ cbt, eliminate = stratum, 
                    family = quasipoisson(), data = data, subset=keep_female)

# Make predictions with all models
pred_male <- crosspred(cbt, model_male, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_female <- crosspred(cbt, model_female, at = seq(minT, maxT, by = 0.1), cen = minT, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Visually compare results of all models
col <- c( "springgreen4","mediumpurple")
par(mfrow = c(1, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_male, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.03)))

lines(pred_female, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.03)))



legend("topleft", c("male","female"), lty=c(1,2,3,4,5,6), lwd=1, col=col,
       bty="n", inset=0.05, y.intersp=1, cex=0.7)

hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.15))


