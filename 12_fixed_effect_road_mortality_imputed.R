source("11_model_setup_imputed.R")

# Model data --------------------------------------------------------------

fit_imputed <- 1:n_imp |>
  map(\(i) {
    all <- fread(paste0(path_data, "final_temp_cluster_subgroup_", i, ".csv"))
    
    # Select all temperature clusters for analysis
    data <- all[cluster_ward_std_6 %in% temp_clusters]
    
    # Specify stratification scheme
    # stratify by city-year-month-dow
    data[, stratum:=factor(paste(salid1, year_month, dow, sep=":"))]
    
    # Keep all strata with nonzero number of road deaths
    data[,  keep:=sum(road)>0, by=stratum]
    
    # Fit model
    model_road <- gnm(road ~ cbt, eliminate = stratum, 
                      family = quasipoisson(), data = data, subset=keep)
    
    # pred_road <- crosspred(cbt, model_road,
    #                        at = seq(minT, maxT, by = .1),
    #                        cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
    
    pred_road_reduced <- crossreduce(cbt, model_road,
                                     at = seq(minT, maxT, by = .1),
                                     cen = medT, bylag = 1, ci.level = 0.95)
    
    # Save coefficients and covariance matrix
    list(coef = coef(model_road),
         vcov = vcov(model_road),
         coef_red = coef(pred_road_reduced),
         vcov_red = vcov(pred_road_reduced))
  }, .progress = TRUE)

# Extract coefficients and covariance matrices
coef_imputed <- map(fit_imputed, \(x) x$coef)
vcov_imputed <- map(fit_imputed, \(x) x$vcov)

coef_imputed_red <- map(fit_imputed, \(x) x$coef_red)
vcov_imputed_red <- map(fit_imputed, \(x) x$vcov_red)

# Pooled results
pooled_results <- apply_rubins_rule(coef_imputed, vcov_imputed)
pooled_results_red <- apply_rubins_rule(coef_imputed_red, vcov_imputed_red)

# Make DLNM predictions with q_bar and T_0
pred_road_bar <- crosspred(cbt,
                           coef = pooled_results$q_bar,
                           vcov = pooled_results$T_0,
                           model.link = "log",
                           at = seq(minT, maxT, by = .1),
                           cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)

pred_road_bar_red <- crossreduce(cbt_red,
                                 coef = pooled_results_red$q_bar,
                                 vcov = pooled_results_red$T_0,
                                 model.link = "log",
                                 at = seq(minT, maxT, by = .1),
                                 cen = medT, bylag = 1, ci.level = 0.95)


# Plot pooled results -----------------------------------------------------

# Plot overal RR and temperature histogram (Pool all coefficients)
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road_bar, "overall", lag=0, cumul=TRUE,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road_bar$predvar<=medT
ind2 <- pred_road_bar$predvar>=medT
lines(pred_road_bar$predvar[ind1],pred_road_bar$allRRfit[ind1],col=4,lwd=1.5)
lines(pred_road_bar$predvar[ind2],pred_road_bar$allRRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60",
     main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))

# Plot overal RR and temperature histogram (Pool reduced coefficients)
par(mfrow = c(2, 1))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_road_bar_red,
     ylim=c(0.8,1.3), ylab="RR", col="brown4", lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha("black", 0.03)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

ind1 <- pred_road_bar_red$predvar<=medT
ind2 <- pred_road_bar_red$predvar>=medT
lines(pred_road_bar_red$predvar[ind1],pred_road_bar_red$RRfit[ind1],col=4,lwd=1.5)
lines(pred_road_bar_red$predvar[ind2],pred_road_bar_red$RRfit[ind2],col=2,lwd=1.5)


hist(data[[Temp_measure]], breaks = 100, probability = TRUE, col = "grey60",
     main = '', xlab = "Temperature (%tile)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.25))
