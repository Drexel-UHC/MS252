################################################################################
# Purpose:
# This script conducts subgroup analysis by age to evaluate
# temperature-related road mortality risks.
#
# It includes:
# - DLNM modeling using the median of 100 imputations.
# - DLNM modeling using imputed datasets (100 imputations).
# - Rubin’s rule to combine imputed model estimates.
# - Plotting male and female pooled exposure-response curves.
#
# Input: Median imputations + 100 imputed age-specific mortality files
# Output: RDS files of pooled coefficients + plots comparing age group-specific effects
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# --- Non-imputed (median) models by age group ----------------------------------
data[,  keep_age1 := sum(deaths_under_9) > 0, by = stratum]
model_age1 <- gnm(deaths_under_9 ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset = keep_age1)

data[,  keep_age2 := sum(deaths_10_19) > 0, by = stratum]
model_age2 <- gnm(deaths_10_19 ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset = keep_age2)

data[,  keep_age3 := sum(deaths_20_34) > 0, by = stratum]
model_age3 <- gnm(deaths_20_34 ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset = keep_age3)

data[,  keep_age4 := sum(deaths_35_64) > 0, by = stratum]
model_age4 <- gnm(deaths_35_64 ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset = keep_age4)

data[,  keep_age5 := sum(deaths_65_plus) > 0, by = stratum]
model_age5 <- gnm(deaths_65_plus ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset = keep_age5)

pred_age1 <- crosspred(cbt, model_age1, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)
pred_age2 <- crosspred(cbt, model_age2, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)
pred_age3 <- crosspred(cbt, model_age3, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)
pred_age4 <- crosspred(cbt, model_age4, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)
pred_age5 <- crosspred(cbt, model_age5, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)

# --- Imputed models: Define analysis function by age group ----------------------
analyze_road_data_age <- function(file_name) {
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path) %>% filter(!grepl("^2020", year))
  
  Temp_measure <- "tmp_pw_percentile"
  minT <- min(data[[Temp_measure]], na.rm = TRUE)
  maxT <- max(data[[Temp_measure]], na.rm = TRUE)
  medT <- median(data[[Temp_measure]], na.rm = TRUE)
  quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
  
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(data[[Temp_measure]], c(10, 75, 90)/100, na.rm = TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group = data$salid1)
  
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  age_cols <- paste0("age", 1:5, "_", sub(".csv", "", file_name))
  
  out <- lapply(age_cols, function(age_col) {
    data[, keep := sum(get(age_col)) > 0, by = stratum]
    model <- gnm(get(age_col) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep)
    pred <- crosspred(cbt, model, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE)
    list(coef = coef(pred), vcov = vcov(pred))
  })
  
  names(out) <- paste0("age", 1:5)
  return(out)
}

# --- Apply function to imputed files and extract results ------------------------
results_age <- lapply(road_files, analyze_road_data_age)
age_groups <- paste0("age", 1:5)

# Extract and pool results using Rubin's rule
for (age in age_groups) {
  coef_list <- lapply(results_age, function(x) x[[age]]$coef)
  vcov_list <- lapply(results_age, function(x) x[[age]]$vcov)
  pooled <- apply_rubin_rule(coef_list, vcov_list)
  saveRDS(pooled$pooled_coef, file = paste0(path_for_pooled, "pooled_", age, "_results_coef.rds"))
  saveRDS(pooled$pooled_vcov, file = paste0(path_for_pooled, "pooled_", age, "_results_vcov.rds"))
}

# --- Reconstruction for plotting ------------------------------------------------
pred_pooled_list <- setNames(lapply(age_groups, function(age) {
  crosspred(cbt, model.link = "log", cum = TRUE, cen = quan01, by = 0.1,
            coef = readRDS(paste0(path_for_pooled, "pooled_", age, "_results_coef.rds")),
            vcov = readRDS(paste0(path_for_pooled, "pooled_", age, "_results_vcov.rds")))
}), age_groups)

# --- Plot results ----------------------------------------------------------------
col <- brewer.pal(n = 21, name = "OrRd")[c(4,5,7,8,9)]

par(mfrow = c(1, 2))
par(mar = c(4, 4, 0.5, 0.5), las = 1, mgp = c(2.5, 1, 0))
plot(pred_pooled_list[['age1']], "overall", lag = 0, cumul = TRUE, ylim = c(0.95, 1.3),
     ylab = "RR", col = col[1], lwd = 1.5, ci = 'bars', lty = 1,
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col[1], 0.01)))
abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")

for (i in 2:5) {
  lines(pred_pooled_list[[age_groups[i]]], "overall", lag = 0, cumul = TRUE, col = col[i],
        lwd = 1.5, lty = 1, ci = "bars", ci.arg = list(col = alpha(col[i], 0.01)))
}

legend("topleft", c("<=9", "10-19", "20-34", "35-64", ">=65"),
       col = col, lty = 1, lwd = 2, bty = "n", cex = 1)

# Overlay median-imputation lines
plot(pred_age1, "overall", lag = 0, cumul = TRUE, ylim = c(0.95, 1.3),
     ylab = "RR", col = col[1], lwd = 1.5, lty = 2, ci = 'bars',
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col[1], 0.01)))
abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")

for (i in 2:5) {
  lines(get(paste0("pred_age", i)), "overall", lag = 0, cumul = TRUE,
        col = col[i], lwd = 1.5, lty = 2, ci = "bars", ci.arg = list(col = alpha(col[i], 0.01)))
}

legend("topleft", c("<=9", "10-19", "20-34", "35-64", ">=65"),
       col = col, lty = 2, lwd = 2, bty = "n", cex = 1)



pred_pooled_list$age5$allRRfit[991]
pred_pooled_list$age5$allRRlow[991]
pred_pooled_list$age5$allRRhigh[991]
