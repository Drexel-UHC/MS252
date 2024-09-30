source("11_model_setup_imputed.R")

# Model data --------------------------------------------------------------

fit_mode_imputed <- 1:n_imp |> 
  map(\(i) {
    all <- fread(paste0(path_data, "final_temp_cluster_subgroup_", i, ".csv"))
    
    # Select all temperature clusters for analysis
    data <- all[cluster_ward_std_6 %in% temp_clusters]
    
    # Specify stratification scheme
    # stratify by city-year-month-dow
    data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]
    
    # Fit model on strata with nonzero number of car deaths
    data[,  keep_veh:=sum(vehicle)>0, by=stratum]
    model_veh <- gnm(vehicle ~ cbt, eliminate = stratum, 
                     family = quasipoisson(), data = data, subset=keep_veh)
    pred_veh_red <- crossreduce(cbt, model_veh,
                                at = seq(minT, maxT, by = .1),
                                cen = medT, bylag = 1, ci.level = 0.95)
    
    # Fit model on strata with nonzero number of motorcycle deaths
    data[,  keep_mot:=sum(motorcycle)>0, by=stratum]
    model_mot <- gnm(motorcycle ~ cbt, eliminate = stratum, 
                     family = quasipoisson(), data = data, subset=keep_mot)
    pred_mot_red <- crossreduce(cbt, model_mot,
                                at = seq(minT, maxT, by = .1),
                                cen = medT, bylag = 1, ci.level = 0.95)
    
    # Fit model on strata with nonzero number of bike deaths
    data[,  keep_bike:=sum(bicycle)>0, by=stratum]
    model_bike <- gnm(bicycle ~ cbt, eliminate = stratum, 
                      family = quasipoisson(), data = data, subset=keep_bike)
    pred_bike_red <- crossreduce(cbt, model_bike,
                                at = seq(minT, maxT, by = .1),
                                cen = medT, bylag = 1, ci.level = 0.95)
    
    # Fit model on strata with nonzero number of pedestrian deaths
    data[,  keep_ped:=sum(ped)>0, by=stratum]
    model_ped <- gnm(ped ~ cbt, eliminate = stratum, 
                     family = quasipoisson(), data = data, subset=keep_ped)
    pred_ped_red <- crossreduce(cbt, model_ped,
                                at = seq(minT, maxT, by = .1),
                                cen = medT, bylag = 1, ci.level = 0.95)
    
    # Save coefficients and covariance matrix
    list(coef_veh = coef(model_veh),
         vcov_veh = vcov(model_veh),
         coef_mot = coef(model_mot),
         vcov_mot = vcov(model_mot),
         coef_bike = coef(model_bike),
         vcov_bike = vcov(model_bike),
         coef_ped = coef(model_ped),
         vcov_ped = vcov(model_ped),
         coef_veh_red = coef(pred_veh_red),
         vcov_veh_red = vcov(pred_veh_red),
         coef_mot_red = coef(pred_mot_red),
         vcov_mot_red = vcov(pred_mot_red),
         coef_bike_red = coef(pred_bike_red),
         vcov_bike_red = vcov(pred_bike_red),
         coef_ped_red = coef(pred_ped_red),
         vcov_ped_red = vcov(pred_ped_red))
  })

# Extract coefficients and covariance matrices
coef_veh <- map(fit_mode_imputed, \(x) x$coef_veh)
vcov_veh <- map(fit_mode_imputed, \(x) x$vcov_veh)
coef_veh_red <- map(fit_mode_imputed, \(x) x$coef_veh_red)
vcov_veh_red <- map(fit_mode_imputed, \(x) x$vcov_veh_red)

coef_mot <- map(fit_mode_imputed, \(x) x$coef_mot)
vcov_mot <- map(fit_mode_imputed, \(x) x$vcov_mot)
coef_mot_red <- map(fit_mode_imputed, \(x) x$coef_mot_red)
vcov_mot_red <- map(fit_mode_imputed, \(x) x$vcov_mot_red)

coef_bike <- map(fit_mode_imputed, \(x) x$coef_bike)
vcov_bike <- map(fit_mode_imputed, \(x) x$vcov_bike)
coef_bike_red <- map(fit_mode_imputed, \(x) x$coef_bike_red)
vcov_bike_red <- map(fit_mode_imputed, \(x) x$vcov_bike_red)

coef_ped <- map(fit_mode_imputed, \(x) x$coef_ped)
vcov_ped <- map(fit_mode_imputed, \(x) x$vcov_ped)
coef_ped_red <- map(fit_mode_imputed, \(x) x$coef_ped_red)
vcov_ped_red <- map(fit_mode_imputed, \(x) x$vcov_ped_red)

# Pooled results
pooled_results_veh <- apply_rubins_rule(coef_veh, vcov_veh)
pooled_results_veh_red <- apply_rubins_rule(coef_veh_red, vcov_veh_red)

pooled_results_mot <- apply_rubins_rule(coef_mot, vcov_mot)
pooled_results_mot_red <- apply_rubins_rule(coef_mot_red, vcov_mot_red)

pooled_results_bike <- apply_rubins_rule(coef_bike, vcov_bike)
pooled_results_bike_red <- apply_rubins_rule(coef_bike_red, vcov_bike_red)

pooled_results_ped <- apply_rubins_rule(coef_ped, vcov_ped)
pooled_results_ped_red <- apply_rubins_rule(coef_ped_red, vcov_ped_red)

# Make predictions with all models
pred_veh <- crosspred(cbt,
                      coef = pooled_results_veh$q_bar,
                      vcov = pooled_results_veh$T_0,
                      model.link = "log",
                      at = seq(minT, maxT, by = .1),
                      cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_veh_red <- crossreduce(cbt_red, 
                            coef = pooled_results_veh_red$q_bar,
                            vcov = pooled_results_veh_red$T_0,
                            model.link = "log",
                            at = seq(minT, maxT, by = .1),
                            cen = medT, bylag = 1, ci.level = 0.95)
pred_mot <- crosspred(cbt,
                      coef = pooled_results_mot$q_bar,
                      vcov = pooled_results_mot$T_0,
                      model.link = "log",
                      at = seq(minT, maxT, by = .1),
                      cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_mot_red <- crossreduce(cbt_red, 
                            coef = pooled_results_mot_red$q_bar,
                            vcov = pooled_results_mot_red$T_0,
                            model.link = "log",
                            at = seq(minT, maxT, by = .1),
                            cen = medT, bylag = 1, ci.level = 0.95)
pred_bike <- crosspred(cbt,
                      coef = pooled_results_bike$q_bar,
                      vcov = pooled_results_bike$T_0,
                      model.link = "log",
                      at = seq(minT, maxT, by = .1),
                      cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_bike_red <- crossreduce(cbt_red, 
                            coef = pooled_results_bike_red$q_bar,
                            vcov = pooled_results_bike_red$T_0,
                            model.link = "log",
                            at = seq(minT, maxT, by = .1),
                            cen = medT, bylag = 1, ci.level = 0.95)
pred_ped <- crosspred(cbt,
                      coef = pooled_results_ped$q_bar,
                      vcov = pooled_results_ped$T_0,
                      model.link = "log",
                      at = seq(minT, maxT, by = .1),
                      cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_ped_red <- crossreduce(cbt_red, 
                            coef = pooled_results_ped_red$q_bar,
                            vcov = pooled_results_ped_red$T_0,
                            model.link = "log",
                            at = seq(minT, maxT, by = .1),
                            cen = medT, bylag = 1, ci.level = 0.95)

# Visually compare results of all models (Pool original coefficients)
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

# Visually compare results of all models (Pool reduced coefficients)
col <- c('firebrick', "darkorange2",'cadetblue','royalblue')
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_veh_red, ylim=c(0.9,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))

lines(pred_mot_red, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))

lines(pred_bike_red, ci="bars",col=col[3], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))

lines(pred_ped_red, ci="bars",col=col[4], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))

legend("topleft", c("car","motorcycle", "bike", "pedestrian"), lty=c(1,2,3,4,5,6), lwd=1, col=col,
       bty="n", inset=0.05, y.intersp=1, cex=0.7)

hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.15))
