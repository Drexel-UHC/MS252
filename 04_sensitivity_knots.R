################################################################################
# Purpose:
# This script performs a sensitivity analysis for the temperature-mortality
# relationship by changing the placement of knots in the exposure-response
# function within a distributed lag non-linear model (DLNM).
#
# Specifically:
# - Uses 100 imputed datasets of daily road deaths.
# - Varies the percentile placement of knots (e.g., 25-75%, 10-90%, etc).
# - For each combination, fits a DLNM and extracts coefficients and variance.
# - Applies Rubin’s rule to combine estimates across imputations.
# - Saves pooled coefficients and variance-covariance matrices for each knot setup.
#
# Input:  100 imputed road death CSVs, temperature data, knot specs.
# Output: Pooled model coefficients (RDS files) for each knot configuration.
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# Load libraries
library(data.table)
library(dlnm)
library(gnm)
library(dplyr)
library(scales)
library(ggplot2)
library(mice)

# Define knot placement specifications (in percentiles)
knot_specs <- list(
  "knot_50" = c(50)
  # "knots_33_67" = c(33.3, 66.7),
  # "knots_25_75" = c(25, 75)
  # "knots_10_90" = c(10, 90),
  # "knots_25_50_75" = c(25, 50, 75)
  # "knots_10_50_90" = c(10, 50, 90)
)

# Base path and file list
base_path <- "/Volumes/TOSHIBA Kai/MS252/Derived Data_20240923_processed_by_imputation/"
road_files <- paste0("road", 1:100, ".csv")

# Function to run analysis for one file and knot configuration
analyze_road_data <- function(file_name, knot_spec, knot_label) {
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path) %>% filter(!grepl("^2020", year))
  
  Temp_measure <- "tmp_pw_percentile"
  minT <- min(data[[Temp_measure]], na.rm = TRUE)
  maxT <- max(data[[Temp_measure]], na.rm = TRUE)
  medT <- median(data[[Temp_measure]], na.rm = TRUE)
  
  knotstmean <- quantile(data[[Temp_measure]], knot_spec / 100, na.rm = TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  
  lagknots <- logknots(3, df = 3)
  cb <- crossbasis(data[[Temp_measure]], lag = 2,
                   argvar = argvartmean,
                   arglag = list(knots = lagknots),
                   group = data$salid1)
  
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  
  model <- gnm(get(road_column) ~ cb, eliminate = stratum,
               family = quasipoisson(), data = data, subset = keep)
  
  pred <- crosspred(cb, model, at = seq(minT, maxT, by = 0.1),
                    cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  return(list(coef = coef(pred), vcov = vcov(pred), knot_label = knot_label))
}

# Run models for each knot setup across all files
results <- list()
for (knot_label in names(knot_specs)) {
  results[[knot_label]] <- lapply(road_files, function(f) {
    analyze_road_data(f, knot_specs[[knot_label]], knot_label)
  }
  )
}

# Rubin's rule function to pool estimates
apply_rubin_rule <- function(coef_list, vcov_list) {
  coef_mat <- do.call(rbind, coef_list)
  vcov_list <- lapply(vcov_list, as.matrix)
  M <- length(vcov_list)
  
  mean_coef <- colMeans(coef_mat)
  var_within <- Reduce("+", vcov_list) / M
  var_between <- Reduce("+", lapply(1:M, function(m) {
    (coef_list[[m]] - mean_coef) %*% t(coef_list[[m]] - mean_coef) / (M - 1)
  })) / M
  total_vcov <- var_within + (1 + 1 / M) * var_between
  
  return(list(pooled_coef = mean_coef, pooled_vcov = total_vcov))
}

# Save pooled results to disk
path_for_pooled <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/sensitivity/"
for (knot_label in names(results)) {
  coef_list <- lapply(results[[knot_label]], `[[`, "coef")
  vcov_list <- lapply(results[[knot_label]], `[[`, "vcov")
  
  pooled_results <- apply_rubin_rule(coef_list, vcov_list)
  
  saveRDS(pooled_results$pooled_coef, paste0(path_for_pooled, knot_label, "_pooled_main_results_coef.rds"))
  saveRDS(pooled_results$pooled_vcov, paste0(path_for_pooled, knot_label, "_pooled_main_results_vcov.rds"))
}


################################################################################
# Purpose:
# The following section visualizes the results of a sensitivity analysis examining how
# different knot placements in the temperature-mortality DLNM affect the
# cumulative exposure-response relationship.
#
# For each knot specification:
# - Loads pooled model coefficients and variances from RDS files.
# - Reconstructs the predicted cumulative relative risk (RR).
# - Plots each curve in a multi-panel figure.
#
# Input:  Pooled coefficient and variance files (from Rubin's rule)
# Output: A multi-panel PNG figure comparing cumulative RR curves across knots
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# Define common components
Temp_measure <- "tmp_pw_percentile"
lagknots <- logknots(3, df = 3)

# Main model (10,75,90) reconstruction
knotstmean <- quantile(data[[Temp_measure]], c(10, 75, 90) / 100, na.rm = TRUE)
argvartmean <- list(fun = "ns", knots = knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, 
                  arglag = list(knots = lagknots), group = data$salid1)
path_main <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/"
main <- crosspred(cbt, model.link = "log",
                  coef = readRDS(paste0(path_main, "pooled_main_results_coef.rds")),
                  vcov = readRDS(paste0(path_main, "pooled_main_results_vcov.rds")),
                  cum = TRUE, cen = quan01, by = 0.1)

# Path to sensitivity results
path_for_pooled <- "/Volumes/TOSHIBA Kai/MS252/Pooled results/sensitivity/"

# Define a list of knot configs and their labels
knot_defs <- list(
  knots_10_50_90 = c(10, 50, 90),
  knots_25_50_75 = c(25, 50, 75),
  knots_10_90 = c(10, 90),
  knots_25_75 = c(25, 75),
  knots_33_67 = c(33.3, 66.7),
  knot_50 = c(50)
)

models <- list()

# Loop over knot definitions to load results and reconstruct predictions
for (k in names(knot_defs)) {
  knotstmean <- quantile(data[[Temp_measure]], knot_defs[[k]] / 100, na.rm = TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, 
                    arglag = list(knots = lagknots), group = data$salid1)
  
  models[[k]] <- crosspred(cbt, model.link = "log",
                           coef = readRDS(paste0(path_for_pooled, k, "_pooled_main_results_coef.rds")),
                           vcov = readRDS(paste0(path_for_pooled, k, "_pooled_main_results_vcov.rds")),
                           cum = TRUE, cen = quan01, by = 0.1)
}

# Define labels for each plot panel
panel_labels <- c(
  "Knots: 10, 50, 90", "Knots: 25, 50, 75", 
  "Knots: 10, 90", "Knots: 25, 75", 
  "Knots: 33, 67", "Knots: 50"
)

# Save to file
png("/Volumes/TOSHIBA Kai/MS252/Figures/figs2.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(4, 2), mar = c(4, 5, 4, 1.5), las = 1, mgp = c(3, 1, 0))

col <- "gray"
  model_names <- names(models)
  
  # Loop through each model and create panel plot
  for (i in seq_along(model_names)) {
    m <- models[[model_names[i]]]
    plot(m, "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4),
         ylab = "RR", col = col, lwd = 1.5, ci = "area", lty = 1,
         xlab = "Temperature (%tile)",
         ci.arg = list(col = alpha(col, 0.3)))
    
    ind2 <- m$predvar >= quan01 & m$predvar <= quan99
    lines(m$predvar[ind2], m$allRRfit[ind2], col = 'firebrick3', lwd = 1.5)
    abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
    mtext(panel_labels[i], side = 3, line = 1, cex = 1.2, font = 2)
  }
  
  # dev.off()
  
  
