### Non-imputed model for age subgroup
data[,  keep_age1:=sum(deaths_under_9)>0, by=stratum]
model_age1 <- gnm(deaths_under_9 ~ cbt, eliminate = stratum, 
               family = quasipoisson(), data = data, subset=keep_age1)

data[,  keep_age2:=sum(deaths_10_19)>0, by=stratum]
model_age2 <- gnm(deaths_10_19 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_age2)

data[,  keep_age3:=sum(deaths_20_34)>0, by=stratum]
model_age3 <- gnm(deaths_20_34 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_age3)

data[,  keep_age4:=sum(deaths_35_64)>0, by=stratum]
model_age4 <- gnm(deaths_35_64 ~ cbt, eliminate = stratum, 
                   family = quasipoisson(), data = data, subset=keep_age4)

data[,  keep_age5:=sum(deaths_65_plus)>0, by=stratum]
model_age5 <- gnm(deaths_65_plus ~ cbt, eliminate = stratum, 
                family = quasipoisson(), data = data, subset=keep_age5)

pred_age1 <- crosspred(cbt, model_age1, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_age2 <- crosspred(cbt, model_age2, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_age3 <- crosspred(cbt, model_age3, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_age4 <- crosspred(cbt, model_age4, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_age5 <- crosspred(cbt, model_age5, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)

# Imputed models by age
# Define function to handle age subgroup analysis for imputed data
analyze_road_data_age <- function(file_name) {
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
  age1_column <- paste0("age1_", sub(".csv", "", file_name))  
  age2_column <- paste0("age2_", sub(".csv", "", file_name))
  age3_column <- paste0("age3_", sub(".csv", "", file_name))  
  age4_column <- paste0("age4_", sub(".csv", "", file_name))
  age5_column <- paste0("age5_", sub(".csv", "", file_name))  
  
  
  ## Age1 deaths analysis
  data[, keep_age1 := sum(get(age1_column)) > 0, by = stratum]
  model_age1 <- gnm(get(age1_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_age1)
  pred_age1 <- crosspred(cbt, model_age1, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_age1 <- coef(pred_age1)
  vcov_age1 <- vcov(pred_age1)
  
  ## Age2 deaths analysis
  data[, keep_age2 := sum(get(age2_column)) > 0, by = stratum]
  model_age2 <- gnm(get(age2_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_age2)
  pred_age2 <- crosspred(cbt, model_age2, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_age2 <- coef(pred_age2)
  vcov_age2 <- vcov(pred_age2)
  
  ## Age3 deaths analysis
  data[, keep_age3 := sum(get(age3_column)) > 0, by = stratum]
  model_age3 <- gnm(get(age3_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_age3)
  pred_age3 <- crosspred(cbt, model_age3, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_age3 <- coef(pred_age3)
  vcov_age3 <- vcov(pred_age3)
  
  ## age4 deaths analysis
  data[, keep_age4 := sum(get(age4_column)) > 0, by = stratum]
  model_age4 <- gnm(get(age4_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_age4)
  pred_age4 <- crosspred(cbt, model_age4, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_age4 <- coef(pred_age4)
  vcov_age4 <- vcov(pred_age4)
  
  ## age5 deaths analysis
  data[, keep_age5 := sum(get(age5_column)) > 0, by = stratum]
  model_age5 <- gnm(get(age5_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_age5)
  pred_age5 <- crosspred(cbt, model_age5, at = seq(minT, maxT, by = 0.1), cen = quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_age5 <- coef(pred_age5)
  vcov_age5 <- vcov(pred_age5)
  
  # Return results for all age groups
  return(list(coef_age1 = coef_age1, vcov_age1 = vcov_age1, 
              coef_age2 = coef_age2, vcov_age2 = vcov_age2,
              coef_age3 = coef_age3, vcov_age3 = vcov_age3, 
              coef_age4 = coef_age4, vcov_age4 = vcov_age4,
              coef_age5 = coef_age5, vcov_age5 = vcov_age5))
}

# Loop through all files for age subgroup analysis
results_age <- lapply(road_files, analyze_road_data_age)

# Define mode groups and result file suffixes
age_groups <- c("age1", "age2", "age3", "age4", "age5")
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"

# Extract coefficients and variance-covariance matrices for all mode groups
coef_list <- setNames(lapply(age_groups, function(age_group) lapply(results_age, function(x) x[[paste0("coef_", age_group)]])), age_groups)
vcov_list <- setNames(lapply(age_groups, function(age_group) lapply(results_age, function(x) x[[paste0("vcov_", age_group)]])), age_groups)

# Apply Rubin's rule and save results for all mode groups
for (age_group in age_groups) {
  # Apply Rubin's rule
  pooled_results <- apply_rubin_rule(coef_list[[age_group]], vcov_list[[age_group]])  # Access by age_group name
  
  # Save RDS files
  saveRDS(pooled_results$pooled_coef, file = paste0(path_for_pooled, "pooled_", age_group, "_results_coef.rds"))
  saveRDS(pooled_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_", age_group, "_results_vcov.rds"))
}



# Reconstruction
pred_pooled_list <- setNames(lapply(age_groups, function(age_group) {
  crosspred(cbt, model.link = "log",
            coef = readRDS(paste0(path_for_pooled, "pooled_", age_group, "_results_coef.rds")),
            vcov = readRDS(paste0(path_for_pooled, "pooled_", age_group, "_results_vcov.rds")),
            cum = TRUE, cen = quan01, by = 0.1)
}), age_groups)

# Plotting
col <- c('darkorange2','darkorange4', "steelblue2","steelblue4",'cadetblue4')
col <- heat.colors(10)[3:8]
col <- brewer.pal(n = 21, name = "OrRd")[c(4,5,7,8,9)]

par(mfrow = c(1, 2))
par(mar=c(4,4,0.5,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_pooled_list[['age1']], "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_pooled_list[['age2']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_pooled_list[['age3']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_pooled_list[['age4']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))
lines(pred_pooled_list[['age5']], "overall", lag=0, cumul=TRUE, ci="bars",col=col[5], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[5], 0.01)))

legend("topleft", c("<=9","10-19","20-34","35-64",">=65"), lty=1, lwd=2, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)

plot(pred_age1, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_age2, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
lines(pred_age3, "overall", lag=0, cumul=TRUE, ci="bars",col=col[3], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[3], 0.01)))
lines(pred_age4, "overall", lag=0, cumul=TRUE, ci="bars",col=col[4], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[4], 0.01)))
lines(pred_age5, "overall", lag=0, cumul=TRUE, ci="bars",col=col[5], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[5], 0.01)))

legend("topleft", c("<=9","10-19","20-34","35-64",">=65"), lty=2, lwd=2, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)



pred_pooled_list[['age1']]$allRRfit[951]
pred_pooled_list[['age1']]$allRRlow[951]
pred_pooled_list[['age1']]$allRRhigh[951]

pred_pooled_list[['age2']]$allRRfit[951]
pred_pooled_list[['age2']]$allRRlow[951]
pred_pooled_list[['age2']]$allRRhigh[951]

pred_pooled_list[['age3']]$allRRfit[951]
pred_pooled_list[['age3']]$allRRlow[951]
pred_pooled_list[['age3']]$allRRhigh[951]

pred_pooled_list[['age4']]$allRRfit[951]
pred_pooled_list[['age4']]$allRRlow[951]
pred_pooled_list[['age4']]$allRRhigh[951]

pred_pooled_list[['age5']]$allRRfit[951]
pred_pooled_list[['age5']]$allRRlow[951]
pred_pooled_list[['age5']]$allRRhigh[951]