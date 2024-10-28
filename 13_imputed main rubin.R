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
base_path <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Derived Data_20240923_processed_by_imputation/"
  

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
  pred_road_reduce <- crossreduce(cbt, model_road, cen=quan01)
  coef_red = coef(pred_road_reduce)
  vcov_red = vcov(pred_road_reduce)
  
  # Return coefficients and standard errors
  return(list(coef = coef, vcov = vcov, coef_red = coef_red, vcov_red = vcov_red))
  
}

analyze_road_data_celcius <- function(file_name) {
  # Read the CSV file
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path)
  
  # Filter data as needed (e.g., covid year)
  data <- data %>% filter(!grepl("^2020", year))
  # Create stratum variable and filter data
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))  # e.g., all_road1 for road1.csv
  # Remove those with no case in strata
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  # Define temperature measure
  Temp_measure_celcius <- "L1ADtemp_pw"
  minT_celcius <- min(data[[Temp_measure_celcius]], na.rm=TRUE)
  maxT_celcius <- max(data[[Temp_measure_celcius]], na.rm=TRUE)
  medT_celcius <- median(data[[Temp_measure_celcius]], na.rm=TRUE)
  quan01_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.01, na.rm = TRUE)
  quan25_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.25, na.rm = TRUE)
  quan75_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.75, na.rm = TRUE)
  quan99_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.99, na.rm = TRUE)
  
  # Create lag knots for crossbasis
  lagknots_celcius <- logknots(3, df = 3)
  knotstmean_celcius <- quantile(data[[Temp_measure_celcius]], c(10,75,90)/100, na.rm=T)
  argvartmean_celcius <- list(fun="ns", knots=knotstmean_celcius)
  cbt_celcius <- crossbasis(data[[Temp_measure_celcius]], lag = 2, argvar = argvartmean_celcius, arglag = list(knots = lagknots_celcius), group=data$salid1)
  # Fit dlnm
  model_road_celcius <- gnm(get(road_column) ~ cbt_celcius, eliminate = stratum, family = quasipoisson(), data = data, subset = keep)
  
  # Predict and generate results
  pred_road_celcius <- crosspred(cbt_celcius, model_road_celcius, at = seq(minT_celcius, maxT_celcius, by = 0.1), cen = medT_celcius, bylag = 1, cumul = TRUE, ci.level = 0.95)
  
  # Extract coefficients
  coef <- coef(pred_road_celcius)
  vcov <- vcov(pred_road_celcius)
  
  # "Reduce" results
  pred_road_reduce_celcius <- crossreduce(cbt_celcius, model_road_celcius, cen=quan01_celcius)
  coef_red = coef(pred_road_reduce_celcius)
  vcov_red = vcov(pred_road_reduce_celcius)
  
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
results_celcius <- lapply(road_files, analyze_road_data_celcius)

# Combine all coefficients and variance-covariance matrices into lists
coef_list <- lapply(results, function(x) x$coef)
vcov_list <- lapply(results, function(x) x$vcov)

coef_list_reduced <- lapply(results, function(x) x$coef_red)
vcov_list_reduced <- lapply(results, function(x) x$vcov_red)

coef_list_celcius <- lapply(results_celcius, function(x) x$coef)
vcov_list_celcius <- lapply(results_celcius, function(x) x$vcov)

coef_list_celcius_reduced <- lapply(results_celcius, function(x) x$coef_red)
vcov_list_celcius_reduced <- lapply(results_celcius, function(x) x$vcov_red)

# Apply Rubin's rule to get pooled coefficients and variance-covariance matrix
pooled_results <- apply_rubin_rule(coef_list, vcov_list)
pooled_results_reduced <- apply_rubin_rule(coef_list_reduced, vcov_list_reduced)
pooled_results_celcius <- apply_rubin_rule(coef_list_celcius, vcov_list_celcius)
pooled_results_celcius_reduced <- apply_rubin_rule(coef_list_celcius_reduced, vcov_list_celcius_reduced)

# Save as RDS
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
saveRDS(pooled_results$pooled_coef, file = paste0(path_for_pooled, "pooled_main_results_coef.rds"))
saveRDS(pooled_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_main_results_vcov.rds"))
saveRDS(pooled_results_reduced$pooled_coef, file = paste0(path_for_pooled, "pooled_main_reduced_results_coef.rds"))
saveRDS(pooled_results_reduced$pooled_vcov, file = paste0(path_for_pooled, "pooled_main_reduced_results_vcov.rds"))
saveRDS(pooled_results_celcius$pooled_coef, file = paste0(path_for_pooled, "pooled_main_celcius_results_coef.rds"))
saveRDS(pooled_results_celcius$pooled_vcov, file = paste0(path_for_pooled, "pooled_main_celcius_results_vcov.rds"))
saveRDS(pooled_results_celcius_reduced$pooled_coef, file = paste0(path_for_pooled, "pooled_main_celcius_reduced_results_coef.rds"))
saveRDS(pooled_results_celcius_reduced$pooled_vcov, file = paste0(path_for_pooled, "pooled_main_celcius_reduced_results_vcov.rds"))


# Reconstruction
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
pred_road_reconstruct_pooled <- crosspred(cbt, model.link="log",
                                          coef = readRDS(paste0(path_for_pooled, "pooled_main_results_coef.rds")),
                                          vcov = readRDS(paste0(path_for_pooled, "pooled_main_results_vcov.rds")),
                                          cum=TRUE,cen=quan01,by=0.1)

pred_road_reconstruct_celcius_pooled <- crosspred(cbt_celcius, model.link="log",
                                                  coef = readRDS(paste0(path_for_pooled, "pooled_main_celcius_results_coef.rds")),
                                                  vcov = readRDS(paste0(path_for_pooled, "pooled_main_celcius_results_vcov.rds")),
                                                  cum=TRUE,cen=quan01_celcius,by=0.1)
# Plotting
col='gray'
par(mfrow = c(2, 2))
par(mar=c(4,5,1,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_road_reconstruct_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.4), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
       xlab="Temperature (%tile) - pooled across 100 imputations", ci.arg=list(col=alpha(col, 0.3)))
#ind1 <- pred_road_reconstruct_pooled$predvar<=medT
ind2 <- pred_road_reconstruct_pooled$predvar >= quan01 & pred_road_reconstruct_pooled$predvar <= quan99
#lines(pred_road_reconstruct_pooled$predvar[ind1],pred_road_reconstruct_pooled$allRRfit[ind1],col='steelblue4',lwd=1.5)
lines(pred_road_reconstruct_pooled$predvar[ind2],pred_road_reconstruct_pooled$allRRfit[ind2],col='firebrick3',lwd=1.5)    
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

plot(pred_road_reconstruct_celcius_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.4), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
       xlab="Temperature (°C) - pooled across 100 imputations", ci.arg=list(col=alpha(col, 0.3)))
#ind1 <- pred_road_reconstruct_celcius_pooled$predvar<=medT_celcius
ind2 <- pred_road_reconstruct_celcius_pooled$predvar >= quan01_celcius & pred_road_reconstruct_celcius_pooled$predvar <= quan99_celcius
#lines(pred_road_reconstruct_celcius_pooled$predvar[ind1],pred_road_reconstruct_celcius_pooled$allRRfit[ind1],col='steelblue4',lwd=1.5)
lines(pred_road_reconstruct_celcius_pooled$predvar[ind2],pred_road_reconstruct_celcius_pooled$allRRfit[ind2],col='firebrick3',lwd=1.5)    
abline(v = c(quan01_celcius,quan25_celcius,medT_celcius,quan75_celcius,quan99_celcius), lty = 4, col = "gray")

plot(pred_road, "overall", lag=0, cumul=TRUE, ylim=c(0.9,1.4), ylab="RR", col=col, lwd=1.5,ci='area',lty=1,
     xlab="Temperature (%tile) - median of 100 imputations", ci.arg=list(col=alpha(col, 0.3)))
#ind1 <- pred_road$predvar>=quan01
ind2 <- pred_road$predvar >= quan01 & pred_road$predvar <= quan99
#lines(pred_road$predvar[ind1],pred_road$allRRfit[ind1],col='steelblue4',lwd=1.5)
lines(pred_road$predvar[ind2],pred_road$allRRfit[ind2],col='firebrick3',lwd=1.5)    
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

hist(data[['L1ADtemp_pw']], breaks = 100, probability = TRUE, col = "grey60", main = '', xlab = "Temperature (°C)", ylab = "Probability Density", yaxt = "n", ylim=c(0,0.15))
abline(v = c(quan01_celcius,quan25_celcius,medT_celcius,quan75_celcius,quan99_celcius), lty = 4, col = "gray")



