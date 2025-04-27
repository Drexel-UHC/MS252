################################################################################
# Purpose:
# This script conducts subgroup analysis by sex (male and female) to evaluate
# temperature-related road mortality risks.
#
# It includes:
# - DLNM modeling using the median of 100 imputations.
# - DLNM modeling using imputed datasets (100 imputations).
# - Rubin’s rule to combine imputed model estimates.
# - Plotting male and female pooled exposure-response curves.
#
# Input: Median imputations + 100 imputed sex-specific mortality files
# Output: RDS files of pooled coefficients + plots comparing male/female effects
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# Median imputations model for sex subgroup (median of 100 imputations)
data[, keep_male := sum(male_deaths) > 0, by = stratum]
model_male <- gnm(male_deaths ~ cbt, eliminate = stratum,
                  family = quasipoisson(), data = data, subset = keep_male)

data[, keep_female := sum(female_deaths) > 0, by = stratum]
model_female <- gnm(female_deaths ~ cbt, eliminate = stratum,
                    family = quasipoisson(), data = data, subset = keep_female)

# Predictions for Median imputations
pred_male <- crosspred(cbt, model_male, at = seq(minT, maxT, by = 0.1),
                       cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_female <- crosspred(cbt, model_female, at = seq(minT, maxT, by = 0.1),
                         cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)

coef_male <- coef(pred_male)
vcov_male <- vcov(pred_male)
coef_female <- coef(pred_female)
vcov_female <- vcov(pred_female)

# Function to perform sex-specific analysis for imputed datasets
analyze_road_data_sex <- function(file_name) {
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path) %>% filter(!grepl("^2020", year))
  
  Temp_measure <- "tmp_pw_percentile"
  minT <- min(data[[Temp_measure]], na.rm = TRUE)
  maxT <- max(data[[Temp_measure]], na.rm = TRUE)
  medT <- median(data[[Temp_measure]], na.rm = TRUE)
  quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
  
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(data[[Temp_measure]], c(10, 75, 90) / 100, na.rm = TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                    argvar = argvartmean, arglag = list(knots = lagknots),
                    group = data$salid1)
  
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  male_column <- paste0("male_", sub(".csv", "", file_name))
  female_column <- paste0("female_", sub(".csv", "", file_name))
  
  data[, keep_male := sum(get(male_column)) > 0, by = stratum]
  model_male <- gnm(get(male_column) ~ cbt, eliminate = stratum,
                    family = quasipoisson(), data = data, subset = keep_male)
  pred_male <- crosspred(cbt, model_male, at = seq(minT, maxT, by = 0.1),
                         cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  data[, keep_female := sum(get(female_column)) > 0, by = stratum]
  model_female <- gnm(get(female_column) ~ cbt, eliminate = stratum,
                      family = quasipoisson(), data = data, subset = keep_female)
  pred_female <- crosspred(cbt, model_female, at = seq(minT, maxT, by = 0.1),
                           cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  return(list(coef_male = coef(pred_male), vcov_male = vcov(pred_male),
              coef_female = coef(pred_female), vcov_female = vcov(pred_female)))
}

# Run across imputations
results_sex <- lapply(road_files, analyze_road_data_sex)

# Separate male/female results
coef_male_list <- lapply(results_sex, `[[`, "coef_male")
vcov_male_list <- lapply(results_sex, `[[`, "vcov_male")
coef_female_list <- lapply(results_sex, `[[`, "coef_female")
vcov_female_list <- lapply(results_sex, `[[`, "vcov_female")

# Apply Rubin’s rule
pooled_male_results <- apply_rubin_rule(coef_male_list, vcov_male_list)
pooled_female_results <- apply_rubin_rule(coef_female_list, vcov_female_list)

# Save pooled estimates
path_for_pooled <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/"
saveRDS(pooled_male_results$pooled_coef, paste0(path_for_pooled, "pooled_male_results_coef.rds"))
saveRDS(pooled_male_results$pooled_vcov, paste0(path_for_pooled, "pooled_male_results_vcov.rds"))
saveRDS(pooled_female_results$pooled_coef, paste0(path_for_pooled, "pooled_female_results_coef.rds"))
saveRDS(pooled_female_results$pooled_vcov, paste0(path_for_pooled, "pooled_female_results_vcov.rds"))

# Reconstruct for plotting
path_for_pooled <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/"
pred_male_pooled <- crosspred(cbt, model.link = "log",
                              coef = readRDS(paste0(path_for_pooled, "pooled_male_results_coef.rds")),
                              vcov = readRDS(paste0(path_for_pooled, "pooled_male_results_vcov.rds")),
                              cum = TRUE, cen = quan01, by = 0.1)

pred_female_pooled <- crosspred(cbt, model.link = "log",
                                coef = readRDS(paste0(path_for_pooled, "pooled_female_results_coef.rds")),
                                vcov = readRDS(paste0(path_for_pooled, "pooled_female_results_vcov.rds")),
                                cum = TRUE, cen = quan01, by = 0.1)

# Plotting
col <- RColorBrewer::brewer.pal(n = 11, name = "Blues")[c(9, 7)]

par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5), las = 1, mgp = c(2.5, 1, 0))

# Plot pooled
plot(pred_male_pooled, "overall", lag = 0, cumul = TRUE, ylim = c(0.95, 1.3),
     ylab = "RR", col = col[1], lwd = 1.5, ci = 'bars', lty = 1,
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col[1], 0.01)))
abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
lines(pred_female_pooled, "overall", lag = 0, cumul = TRUE, col = col[2],
      lty = 1, lwd = 1.5, ci = "bars", ci.arg = list(col = alpha(col[2], 0.01)))
legend("topleft", legend = c("male", "female"), col = col, lty = 1, lwd = 1,
       bty = "n", inset = 0, y.intersp = 0.8, cex = 1.0)

# Plot Median imputations
plot(pred_male, "overall", lag = 0, cumul = TRUE, ylim = c(0.95, 1.3),
     ylab = "RR", col = col[1], lwd = 1.5, ci = 'bars', lty = 2,
     xlab = "Temperature (%tile)", ci.arg = list(col = alpha(col[1], 0.01)))
abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
lines(pred_female, "overall", lag = 0, cumul = TRUE, col = col[2],
      lty = 2, lwd = 1.5, ci = "bars", ci.arg = list(col = alpha(col[2], 0.01)))
legend("topleft", legend = c("male", "female"), col = col, lty = 2, lwd = 1,
       bty = "n", inset = 0, y.intersp = 0.8, cex = 1.0)

# Inspect RR values
pred_male_pooled$allRRfit[991]
pred_male_pooled$allRRlow[991]
pred_male_pooled$allRRhigh[991]
pred_female_pooled$allRRfit[991]
pred_female_pooled$allRRlow[991]
pred_female_pooled$allRRhigh[991]
