library(data.table)
library(dlnm)
library(gnm)
library(dplyr)
library(scales) # For the alpha function
library(ggplot2)
library(mice)

# check if road1-road100 have same dimensions and have 272 salid1


# Define the range of files to process
road_files <- paste0("road", 1:100, ".csv")  # Adjust this as needed for your specific files

# Define the base directory path
base_path <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS 252_all_imputed/Derived Data_processed_by_imputation/"

# Function to perform analysis on each file
analyze_road_data <- function(file_name) {
  # Read the CSV file
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path)
  
  # Filter data as needed (e.g., covid year)
  data <- data %>% filter(!grepl("^2020", year))
  
  # Define temperature measure
  Temp_measure <- "tmp_pw_percentile"
  #Temp_measure <- "L1ADtemp_pw"
  
  # Calculate quantiles and other statistics
  minT <- min(data[[Temp_measure]], na.rm=TRUE)
  maxT <- max(data[[Temp_measure]], na.rm=TRUE)
  medT <- median(data[[Temp_measure]], na.rm=TRUE)
  quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
  quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
  quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
  quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)
  
  # Create lag knots for crossbasis
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  
  # Create crossbasis object
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
  
  # Create stratum variable and filter data
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))  # e.g., all_road1 for road1.csv
  
  # Remove those with no case in strata
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  
  # Fit dlnm
  model_road <- gnm(get(road_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep)
  
  # Predict and generate results
  pred_road <- crosspred(cbt, model_road, at = seq(minT, maxT, by = 0.1), cen = medT, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  # Extract coefficients
  coef <- coef(pred_road)
  vcov <- vcov(pred_road)
  
  # "Reduce" results
  pred_road_reduce <- crossreduce(cbt, model_road)
  coef_red = coef(pred_road_reduce)
  vcov_red = vcov(pred_road_reduce)
  
  # Return coefficients and standard errors
  return(list(coef = coef, vcov = vcov, coef_red = coef_red, vcov_red = vcov_red))
  
}

# Function to apply Rubin's rule and get pooled estimates and variance-covariance matrix
apply_rubin_rule <- function(coef_list, vcov_list) {
  # Convert lists to matrices
  coef_mat <- do.call(rbind, coef_list)
  vcov_list <- lapply(vcov_list, function(x) as.matrix(x))  # Ensure vcov_list is a list of matrices
  
  # Number of imputations
  M <- length(vcov_list)
  
  # Calculate the mean coefficients
  mean_coef <- colMeans(coef_mat)
  
  # Calculate the within-imputation variance-covariance matrix
  var_within <- Reduce("+", vcov_list) / M
  
  # Calculate the between-imputation variance-covariance matrix
  var_between <- Reduce("+", lapply(1:M, function(m) {
    (coef_list[[m]] - mean_coef) %*% t(coef_list[[m]] - mean_coef) / (M - 1)
  })) / M
  
  # Calculate the pooled variance-covariance matrix
  total_vcov <- var_within + (1 + 1/M) * var_between
  
  # Return pooled coefficients and variance-covariance matrix
  return(list(
    pooled_coef = mean_coef,
    pooled_vcov = total_vcov
  ))
}

# Loop through all files, run the analysis, and store coefficients and variance-covariance matrices in lists
results <- lapply(road_files, analyze_road_data)

# Combine all coefficients and variance-covariance matrices into lists
coef_list <- lapply(results, function(x) x$coef)
vcov_list <- lapply(results, function(x) x$vcov)

coef_list_reduced <- lapply(results, function(x) x$coef_red)
vcov_list_reduced <- lapply(results, function(x) x$vcov_red)

# Apply Rubin's rule to get pooled coefficients and variance-covariance matrix
pooled_results <- apply_rubin_rule(coef_list, vcov_list)
pooled_results_reduced <- apply_rubin_rule(coef_list_reduced, vcov_list_reduced)

# Print or save the pooled estimates and variance-covariance matrix
print(pooled_results$pooled_coef)
print(pooled_results$pooled_vcov)
print(pooled_results_reduced$pooled_coef)
print(pooled_results_reduced$pooled_vcov)

