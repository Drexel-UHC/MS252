### Non-imputed model for mode subgroup
data[,  keep_vehicle:=sum(median_vehicle_round)>0, by=stratum]
model_vehicle <- gnm(median_vehicle_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_vehicle)

data[,  keep_motorcycle:=sum(median_motorcycle_round)>0, by=stratum]
model_motorcycle <- gnm(median_motorcycle_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_motorcycle)

data[,  keep_bicycle:=sum(median_bicycle_round)>0, by=stratum]
model_bicycle <- gnm(median_bicycle_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_bicycle)

data[,  keep_ped:=sum(median_ped_round)>0, by=stratum]
model_ped <- gnm(median_ped_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_ped)


pred_vehicle <- crosspred(cbt, model_vehicle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_motorcycle <- crosspred(cbt, model_motorcycle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_bicycle <- crosspred(cbt, model_bicycle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_ped <- crosspred(cbt, model_ped, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Imputed models by mode
# Define function to handle mode subgroup analysis for imputed data
analyze_road_data_mode <- function(file_name) {
  # Read the CSV file
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path)
  
  # Filter data (exclude covid year)
  data <- data %>% filter(!grepl("^2020", year))
  
  # Define temperature measure
  Temp_measure <- "tmp_pw_percentile"
  
  # Calculate quantiles and other statistics
  minT <- min(data[[Temp_measure]], na.rm=TRUE)
  maxT <- max(data[[Temp_measure]], na.rm=TRUE)
  medT <- median(data[[Temp_measure]], na.rm=TRUE)
  quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
  quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
  quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
  quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)
  
  # Create crossbasis object
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(data[[Temp_measure]], c(10, 75, 90)/100, na.rm=TRUE)
  argvartmean <- list(fun = "ns", knots = knotstmean)
  cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
  
  # Stratum variable
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  
  # Adjust for specific death columns based on file name
  vehicle_column <- paste0("vehicle_", sub(".csv", "", file_name))  
  motorcycle_column <- paste0("motorcycle_", sub(".csv", "", file_name))
  bicycle_column <- paste0("bicycle_", sub(".csv", "", file_name))  
  ped_column <- paste0("ped_", sub(".csv", "", file_name))

  
  ## vehicle deaths analysis
  data[, keep_vehicle := sum(get(vehicle_column)) > 0, by = stratum]
  model_vehicle <- gnm(get(vehicle_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_vehicle)
  pred_vehicle <- crosspred(cbt, model_vehicle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_vehicle <- coef(pred_vehicle)
  vcov_vehicle <- vcov(pred_vehicle)
  
  ## motorcycle deaths analysis
  data[, keep_motorcycle := sum(get(motorcycle_column)) > 0, by = stratum]
  model_motorcycle <- gnm(get(motorcycle_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_motorcycle)
  pred_motorcycle <- crosspred(cbt, model_motorcycle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_motorcycle <- coef(pred_motorcycle)
  vcov_motorcycle <- vcov(pred_motorcycle)
  
  ## bicycle deaths analysis
  data[, keep_bicycle := sum(get(bicycle_column)) > 0, by = stratum]
  model_bicycle <- gnm(get(bicycle_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_bicycle)
  pred_bicycle <- crosspred(cbt, model_bicycle, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_bicycle <- coef(pred_bicycle)
  vcov_bicycle <- vcov(pred_bicycle)
  
  ## ped deaths analysis
  data[, keep_ped := sum(get(ped_column)) > 0, by = stratum]
  model_ped <- gnm(get(ped_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_ped)
  pred_ped <- crosspred(cbt, model_ped, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_ped <- coef(pred_ped)
  vcov_ped <- vcov(pred_ped)
  
  # Return results for all mode groups
  return(list(coef_vehicle = coef_vehicle, vcov_vehicle = vcov_vehicle, 
              coef_motorcycle = coef_motorcycle, vcov_motorcycle = vcov_motorcycle,
              coef_bicycle = coef_bicycle, vcov_bicycle = vcov_bicycle, 
              coef_ped = coef_ped, vcov_ped = vcov_ped))
}

# Loop through all files for mode subgroup analysis
results_mode <- lapply(road_files, analyze_road_data_mode)

# Define mode groups and result file suffixes
modes <- c("vehicle", "motorcycle", "bicycle", "ped")
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
# Extract coefficients and variance-covariance matrices for all mode groups
coef_list <- setNames(lapply(modes, function(mode) lapply(results_mode, function(x) x[[paste0("coef_", mode)]])), modes)
vcov_list <- setNames(lapply(modes, function(mode) lapply(results_mode, function(x) x[[paste0("vcov_", mode)]])), modes)

# Apply Rubin's rule and save results for all mode groups
for (mode in modes) {
  # Apply Rubin's rule
  pooled_results <- apply_rubin_rule(coef_list[[mode]], vcov_list[[mode]])  # Access by mode name
  
  # Save RDS files
  saveRDS(pooled_results$pooled_coef, file = paste0(path_for_pooled, "pooled_", mode, "_results_coef.rds"))
  saveRDS(pooled_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_", mode, "_results_vcov.rds"))
}


# Reconstruction
modes <- c("vehicle", "motorcycle", "bicycle", "ped")
pred_pooled_list <- setNames(lapply(modes, function(mode) {
  crosspred(cbt, model.link = "log",
            coef = readRDS(paste0(path_for_pooled, "pooled_", mode, "_results_coef.rds")),
            vcov = readRDS(paste0(path_for_pooled, "pooled_", mode, "_results_vcov.rds")),
            cum = TRUE, cen = quan01, by = 0.1)
}), modes)

# Plotting
col <- c("firebrick2",'darkorange2','cadetblue4', 'steelblue3')
col <- brewer.pal(n = 11, name = "RdYlBu")[c(1,3,9,11)]
par(mfrow = c(1, 2))
par(mar=c(4,4,0.5,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_pooled_list[['motorcycle']], "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_pooled_list[['bicycle']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_pooled_list[['vehicle']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_pooled_list[['ped']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))
legend("topleft", c("motorcycle","bicycle","vehicle","ped"), lty=1, lwd=1, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)

plot(pred_motorcycle, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_bicycle, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_vehicle, "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_ped, "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))

legend("topleft", c("motorcycle","bicycle","vehicle","ped"), lty=2, lwd=1, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)




pred_pooled_list[['motorcycle']]$allRRfit[951]
pred_pooled_list[['motorcycle']]$allRRlow[951]
pred_pooled_list[['motorcycle']]$allRRhigh[951]

pred_pooled_list[['bicycle']]$allRRfit[951]
pred_pooled_list[['bicycle']]$allRRlow[951]
pred_pooled_list[['bicycle']]$allRRhigh[951]


pred_pooled_list[['vehicle']]$allRRfit[951]
pred_pooled_list[['vehicle']]$allRRlow[951]
pred_pooled_list[['vehicle']]$allRRhigh[951]

pred_pooled_list[['ped']]$allRRfit[951]
pred_pooled_list[['ped']]$allRRlow[951]
pred_pooled_list[['ped']]$allRRhigh[951]