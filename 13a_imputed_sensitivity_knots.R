library(data.table)
library(dlnm)
library(gnm)
library(dplyr)
library(scales)
library(ggplot2)
library(mice)

# Define knot specifications
knot_specs <- list(
  "knot_50" = c(50),
  "knots_33_67" = c(33.3, 66.7),
  "knots_25_75" = c(25, 75),
  "knots_10_90" = c(10, 90),
  "knots_25_50_75" = c(25, 50, 75),
  "knots_10_50_90" = c(10, 50, 90)
  # "knots_10_75_90" = c(10, 75, 90) #main effects
)

# Define the base directory path
base_path <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Derived Data_20240923_processed_by_imputation/"

# Define road files
road_files <- paste0("road", 1:100, ".csv")

# Function to perform analysis
analyze_road_data <- function(file_name, knot_spec, knot_label) {
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path)
  
  # Filter out COVID year
  data <- data %>% filter(!grepl("^2020", year))
  
  Temp_measure <- "tmp_pw_percentile"
  
  # Compute relevant quantiles
  minT <- min(data[[Temp_measure]], na.rm=TRUE)
  maxT <- max(data[[Temp_measure]], na.rm=TRUE)
  medT <- median(data[[Temp_measure]], na.rm=TRUE)
  
  # Compute knots
  knotstmean <- quantile(data[[Temp_measure]], knot_spec / 100, na.rm=TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  
  # Create crossbasis
  lagknots <- logknots(3, df = 3)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
  
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  
  # Fit DLNM model
  model_road <- gnm(get(road_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep)
  
  # Predict results
  pred_road <- crosspred(cbt, model_road, at = seq(minT, maxT, by = 0.1), cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  coef <- coef(pred_road)
  vcov <- vcov(pred_road)
  
  return(list(coef = coef, vcov = vcov, knot_label = knot_label))
}

# Run analyses for all files and knot specifications
results <- list()
for (knot_label in names(knot_specs)) {
  results[[knot_label]] <- lapply(road_files, function(f) analyze_road_data(f, knot_specs[[knot_label]], knot_label))
}

# Apply Rubin's rule for pooled estimates
apply_rubin_rule <- function(coef_list, vcov_list) {
  coef_mat <- do.call(rbind, coef_list)
  vcov_list <- lapply(vcov_list, function(x) as.matrix(x))
  M <- length(vcov_list)
  
  mean_coef <- colMeans(coef_mat)
  var_within <- Reduce("+", vcov_list) / M
  var_between <- Reduce("+", lapply(1:M, function(m) {
    (coef_list[[m]] - mean_coef) %*% t(coef_list[[m]] - mean_coef) / (M - 1)
  })) / M
  total_vcov <- var_within + (1 + 1/M) * var_between
  
  return(list(pooled_coef = mean_coef, pooled_vcov = total_vcov))
}

# Store results
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/sensitivity/"
for (knot_label in names(results)) {
  coef_list <- lapply(results[[knot_label]], function(x) x$coef)
  vcov_list <- lapply(results[[knot_label]], function(x) x$vcov)
  
  pooled_results <- apply_rubin_rule(coef_list, vcov_list)
  saveRDS(pooled_results$pooled_coef, file = paste0(path_for_pooled, knot_label, "_pooled_main_results_coef.rds"))
  saveRDS(pooled_results$pooled_vcov, file = paste0(path_for_pooled, knot_label, "_pooled_main_results_vcov.rds"))
}
