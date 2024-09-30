### Non-imputed model for sex subgroup

data[,  keep_male:=sum(male_deaths)>0, by=stratum]
model_male <- gnm(male_deaths ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep_male)

data[,  keep_female:=sum(female_deaths)>0, by=stratum]
model_female <- gnm(female_deaths ~ cbt, eliminate = stratum, 
                    family = quasipoisson(), data = data, subset=keep_female)

pred_male <- crosspred(cbt, model_male, at = seq(minT, maxT, by = 0.1), cen=quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
pred_female <- crosspred(cbt, model_female, at = seq(minT, maxT, by = 0.1), cen=quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)

coef_male = coef(pred_male)
vcov_male = vcov(pred_male)
coef_female = coef(pred_female)
vcov_female = vcov(pred_female)

# Imputed models by sex
# Define function to handle sex subgroup analysis for imputed data
analyze_road_data_sex <- function(file_name) {
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
  
  # Adjust for male and female specific death columns based on file name
  male_column <- paste0("male_", sub(".csv", "", file_name))  # e.g., male_road1 for road1.csv
  female_column <- paste0("female_", sub(".csv", "", file_name))  # e.g., female_road1 for road1.csv
  
  ## Male deaths analysis
  data[, keep_male := sum(get(male_column)) > 0, by = stratum]
  model_male <- gnm(get(male_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_male)
  pred_male <- crosspred(cbt, model_male, at = seq(minT, maxT, by = 0.1), cen=quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_male <- coef(pred_male)
  vcov_male <- vcov(pred_male)
  
  ## Female deaths analysis
  data[, keep_female := sum(get(female_column)) > 0, by = stratum]
  model_female <- gnm(get(female_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep_female)
  pred_female <- crosspred(cbt, model_female, at = seq(minT, maxT, by = 0.1), cen=quan01, bylag = 1, cumul = TRUE, ci.level = 0.95)
  coef_female <- coef(pred_female)
  vcov_female <- vcov(pred_female)
  
  # Return results for both male and female analyses
  return(list(coef_male = coef_male, vcov_male = vcov_male, coef_female = coef_female, vcov_female = vcov_female))
}

# Loop through all files for both male and female subgroup analysis
results_sex <- lapply(road_files, analyze_road_data_sex)

# Extract coefficients and variance-covariance matrices for males and females
coef_male_list <- lapply(results_sex, function(x) x$coef_male)
vcov_male_list <- lapply(results_sex, function(x) x$vcov_male)
coef_female_list <- lapply(results_sex, function(x) x$coef_female)
vcov_female_list <- lapply(results_sex, function(x) x$vcov_female)

# Apply Rubin's rule to get pooled estimates for males and females
pooled_male_results <- apply_rubin_rule(coef_male_list, vcov_male_list)
pooled_female_results <- apply_rubin_rule(coef_female_list, vcov_female_list)

# Save as RDS
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
saveRDS(pooled_male_results$pooled_coef, file = paste0(path_for_pooled, "pooled_male_results_coef.rds"))
saveRDS(pooled_male_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_male_results_vcov.rds"))
saveRDS(pooled_female_results$pooled_coef, file = paste0(path_for_pooled, "pooled_female_results_coef.rds"))
saveRDS(pooled_female_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_female_results_vcov.rds"))

# Reconstruction
pred_male_pooled <- crosspred(cbt, model.link="log",
                              coef = readRDS(paste0(path_for_pooled, "pooled_male_results_coef.rds")),
                              vcov = readRDS(paste0(path_for_pooled, "pooled_male_results_vcov.rds")),
                              cum=TRUE, cen=quan01, by=0.1)

pred_female_pooled <- crosspred(cbt, model.link="log",
                                coef = readRDS(paste0(path_for_pooled, "pooled_female_results_coef.rds")),
                                vcov = readRDS(paste0(path_for_pooled, "pooled_female_results_vcov.rds")),
                                cum=TRUE, cen=quan01, by=0.1)

# Plotting
col <- c( "darkorange2","steelblue3")
col <- brewer.pal(n = 11, name = "Blues")[c(9,7)]

par(mfrow = c(1, 2))
par(mar=c(4,4,0.5,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_male_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_female_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))

legend("topleft", c("male","female"), lty=1, lwd=1, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)

plot(pred_male, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")

lines(pred_female, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))

legend("topleft", c("male","female"), lty=2, lwd=1, col=col,
       bty="n", inset=0, y.intersp=0.8, cex=1.0)


pred_male_pooled$allRRfit[951]
pred_male_pooled$allRRlow[951]
pred_male_pooled$allRRhigh[951]

pred_female_pooled$allRRfit[951]
pred_female_pooled$allRRlow[951]
pred_female_pooled$allRRhigh[951]