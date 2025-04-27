################################################################################
# Purpose:
# This script performs interaction analyses between temperature and multiple
# effect modifiers (e.g., temperature mean, SD, traffic time, street length)
# using imputed datasets for road-traffic deaths. For each modifier:
# - Fits DLNM model with interaction terms at low (10th percentile) and high (90th)
# - Applies Rubin's rule to pool estimates across 100 imputations
# - Outputs pooled coefficients and VCOVs as RDS
# - Reconstructs predictions for plotting
# - Use median-imputation model for test of effect modification
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# ---- Setup ----
path_for_pooled <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/"
modifiers <- c("mean", "std", "BECADSTTLGAVGL1AD", "BECURBAVGTRAFTIMEL1AD")

# ---- Function: Analyze one imputed file with one modifier ----
analyze_road_data_int <- function(file_name, modifier) {
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path) %>% filter(!grepl("^2020", year))
  
  Temp_measure <- "tmp_pw_percentile"
  minT <- min(data[[Temp_measure]], na.rm=TRUE)
  maxT <- max(data[[Temp_measure]], na.rm=TRUE)
  medT <- median(data[[Temp_measure]], na.rm=TRUE)
  quan01 <- quantile(data[[Temp_measure]], 0.01, na.rm=TRUE)
  quan25 <- quantile(data[[Temp_measure]], 0.25, na.rm=TRUE)
  quan75 <- quantile(data[[Temp_measure]], 0.75, na.rm=TRUE)
  quan99 <- quantile(data[[Temp_measure]], 0.99, na.rm=TRUE)
  
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean,
                    arglag = list(knots = lagknots), group = data$salid1)
  
  # Interaction terms
  intval <- quantile(data[[modifier]], c(0.1, 0.9))
  cbint1 <- cbt * (data[[modifier]] - intval[1])
  cbint2 <- cbt * (data[[modifier]] - intval[2])
  
  # Model fitting
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  
  model_road <- gnm(get(road_column) ~ cbt, eliminate = stratum,
                    family = quasipoisson(), data = data, subset = keep)
  modint1 <- update(model_road, . ~ . + cbint1)
  modint2 <- update(model_road, . ~ . + cbint2)
  
  # Predictions
  pred_road_int1 <- crosspred(cbt, modint1, cen=quan01)
  pred_road_int2 <- crosspred(cbt, modint2, cen=quan01)
  
  return(list(coef_int1 = coef(pred_road_int1), vcov_int1 = vcov(pred_road_int1),
              coef_int2 = coef(pred_road_int2), vcov_int2 = vcov(pred_road_int2)))
}

# ---- Loop: Apply to all modifiers and all imputed datasets ----
all_results <- list()
for (modifier in modifiers) {
  results <- lapply(road_files, function(f) analyze_road_data_int(f, modifier))
  
  coef_int1_list <- lapply(results, `[[`, "coef_int1")
  vcov_int1_list <- lapply(results, `[[`, "vcov_int1")
  coef_int2_list <- lapply(results, `[[`, "coef_int2")
  vcov_int2_list <- lapply(results, `[[`, "vcov_int2")
  
  pooled_int1 <- apply_rubin_rule(coef_int1_list, vcov_int1_list)
  pooled_int2 <- apply_rubin_rule(coef_int2_list, vcov_int2_list)
  
  saveRDS(pooled_int1$pooled_coef, paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
  saveRDS(pooled_int1$pooled_vcov, paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
  saveRDS(pooled_int2$pooled_coef, paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
  saveRDS(pooled_int2$pooled_vcov, paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))
  
  all_results[[modifier]] <- list(int1 = pooled_int1, int2 = pooled_int2)
}

# ---- Reconstruction for plotting ----
modifier <- "mean"  # example modifier to plot
coef_int1 <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1 <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2 <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2 <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef=coef_int1, vcov=vcov_int1,
                              cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log", coef=coef_int2, vcov=vcov_int2,
                              cum=TRUE, cen=quan01, by=0.1)

# ---- Plotting ----
col <- brewer.pal(n = 9, name = "Set1")[c(2,1)]
par(mfrow = c(1, 2))
par(mar=c(4,3.5,0.5,0.5), las=1, mgp=c(2.5,1,0))

plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1],
     lwd=1.5, ci='bars', lty=1, xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, col=col[2], lwd=1.5, lty=1, ci='bars',
      ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty=4, col="gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col, bty="n", inset=0.01, y.intersp=0.9, cex=1.0)



# ---- Test of effect modifcation ----
data[, year := as.numeric(format(as.Date(allDate), "%Y"))]
modifier <- 'BECURBAVGTRAFTIMEL1AD' #change to other modifiers in c("mean", "std", "BECADSTTLGAVGL1AD", "BECURBAVGTRAFTIMEL1AD")
intval <- quantile(data[[modifier]], c(0.10, 0.90))
cbint1 <- cbt * (data[[modifier]] - intval[1])
cbint2 <- cbt * (data[[modifier]] - intval[2])
modint1 <- update(model_road, .~. + cbint1)
modint2 <- update(model_road, .~. + cbint2)
anova(model_road, modint1, test="Chisq")
pred_int1 <- crosspred(cbt, modint1, cen=quan01, cum=TRUE, by=0.1)
pred_int2 <- crosspred(cbt, modint2, cen=quan01, cum=TRUE, by=0.1)
