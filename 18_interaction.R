# imputed models
# Function to perform analysis on each file
analyze_road_data_int <- function(file_name, modifier) {
  # Read the CSV file
  modifier = modifier
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
  intval <- quantile(data[[modifier]], c(0.1, 0.9))
  cbint1 <- cbt * (data[[modifier]] - intval[1])
  cbint2 <- cbt * (data[[modifier]] - intval[2])
  
  # Create stratum variable and filter data
  data[, stratum := factor(paste(salid1, year_month_dow, sep = ":"))]
  road_column <- paste0("all_", sub(".csv", "", file_name))  # e.g., all_road1 for road1.csv
  
  # Remove those with no case in strata
  data[, keep := sum(get(road_column)) > 0, by = stratum]
  
  # Fit dlnm
  model_road <- gnm(get(road_column) ~ cbt, eliminate = stratum, family = quasipoisson(), data = data, subset = keep)
  
  modint1 <- update(model_road, .~. + cbint1)
  modint2 <- update(model_road, .~. + cbint2)
  
  # Predict and generate results
  pred_road_int1 <- crosspred(cbt, modint1, cen=quan01)
  pred_road_int2 <- crosspred(cbt, modint2, cen=quan01)
  
  # Extract coefficients
  coef_int1 <- coef(pred_road_int1)
  vcov_int1 <- vcov(pred_road_int1)
  coef_int2 <- coef(pred_road_int2)
  vcov_int2 <- vcov(pred_road_int2)
  
  # Return coefficients and standard errors
  return(list(coef_int1 = coef_int1, vcov_int1 = vcov_int1, 
              coef_int2 = coef_int2, vcov_int2 = vcov_int2))
  
}
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"

# List of modifiers
modifiers <- c("mean", "std", "BECADSTTLGAVGL1AD", "BECURBAVGTRAFTIMEL1AD") #BECADINTDENSL1AD intersection; #'BECADINTDENS3L1AD' #3-leg; 'BECADINTDENS4L1AD' #4-leg; 'BECADSTTPNODEAVGL1AD' #street per node; 'BECADCRCTYAVGL1AD' #Circuity average; 'BECURBTRVDELAYTIMEL1AD'

# Prepare an empty list to store results for all analyses
all_results <- list()

# Loop through each modifier and perform the analysis
for (modifier in modifiers) {
  # Perform analysis on each file and store results in a list
  results <- lapply(road_files, function(file_name) analyze_road_data_int(file_name, modifier))
  
  # Extract coefficients and variance-covariance matrices
  coef_int1_list <- lapply(results, function(x) x$coef_int1)
  vcov_int1_list <- lapply(results, function(x) x$vcov_int1)
  coef_int2_list <- lapply(results, function(x) x$coef_int2)
  vcov_int2_list <- lapply(results, function(x) x$vcov_int2)
  
  # Apply Rubin's rule
  pooled_int1_results <- apply_rubin_rule(coef_int1_list, vcov_int1_list)
  pooled_int2_results <- apply_rubin_rule(coef_int2_list, vcov_int2_list)
  
  # Save results to RDS files
  saveRDS(pooled_int1_results$pooled_coef, file = paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
  saveRDS(pooled_int1_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
  saveRDS(pooled_int2_results$pooled_coef, file = paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
  saveRDS(pooled_int2_results$pooled_vcov, file = paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))
  
  # Store the results in the main list
  all_results[[modifier]] <- list(int1 = pooled_int1_results, int2 = pooled_int2_results)
}

# Reconstruction by modifier
modifier = 'BECURBAVGTRAFTIMEL1AD' #mean #std #BECADSTTLGAVGL1AD #BECURBAVGTRAFTIMEL1AD
coef_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_coef.rds"))
vcov_int1_pooled <- readRDS(paste0(path_for_pooled, "pooled_int1_", modifier, "_results_vcov.rds"))
coef_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_coef.rds"))
vcov_int2_pooled <- readRDS(paste0(path_for_pooled, "pooled_int2_", modifier, "_results_vcov.rds"))

pred_int1_pooled <- crosspred(cbt, model.link="log", coef = coef_int1_pooled,vcov = vcov_int1_pooled, cum=TRUE, cen=quan01, by=0.1)
pred_int2_pooled <- crosspred(cbt, model.link="log",coef = coef_int2_pooled, vcov = vcov_int2_pooled, cum=TRUE, cen=quan01, by=0.1)

### unimputed model for mode subgroup
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


#Plotting
col <- c( 'royalblue', 'tomato')
col <- c("steelblue4", "firebrick3")
par(mfrow = c(1, 2))
par(mar=c(4,3.5,0.5,0.5), las=1, mgp=c(2.5,1,0))
plot(pred_int1_pooled, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=1,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
lines(pred_int2_pooled, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=1,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
legend("topleft", c("lower","higher"), lty=1, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)

plot(pred_int1, "overall", lag=0, cumul=TRUE, ylim=c(0.95,1.3), ylab="RR", col=col[1], lwd=1.5,ci='bars',lty=2,
     xlab="Temperature (%tile)", ci.arg=list(col=alpha(col[1], 0.01)))
abline(v = c(quan01,quan25,medT,quan75,quan99), lty = 4, col = "gray")
lines(pred_int2, "overall", lag=0, cumul=TRUE, ci="bars",col=col[2], lty=2,
      lwd=1.5, ci.arg=list(col=alpha(col[2], 0.01)))
legend("topleft", c("lower","higher"), lty=2, lwd=1, col=col,
       bty="n", inset=0.01, y.intersp=0.9, cex=1.0)
