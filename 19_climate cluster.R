# Function to perform analysis on each file
analyze_road_data_cluster <- function(file_name, cluster) {
  # Read the CSV file
  cluster = cluster
  
  file_path <- paste0(base_path, file_name)
  data <- fread(file_path)
  
  data <- data[data[['cluster_ward_std_6']] %in% c(cluster)] #####################################
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

# Define the number of clusters
num_clusters <- 6

# Loop through each cluster to run the analysis and store results
results_clusters <- vector("list", num_clusters)
coef_lists <- vector("list", num_clusters)
vcov_lists <- vector("list", num_clusters)
pooled_results <- vector("list", num_clusters)

for (i in 1:num_clusters) {
  # Run the analysis for each cluster
  results_clusters[[i]] <- lapply(road_files, function(file_name) analyze_road_data_cluster(file_name, i))
  
  # Store coefficients and variance-covariance matrices
  coef_lists[[i]] <- lapply(results_clusters[[i]], function(x) x$coef)
  vcov_lists[[i]] <- lapply(results_clusters[[i]], function(x) x$vcov)
  
  # Apply Rubin's rule
  pooled_results[[i]] <- apply_rubin_rule(coef_lists[[i]], vcov_lists[[i]])
}

# Save the pooled coefficients to RDS files
path_for_pooled <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Pooled results/"
for (i in 1:num_clusters) {
  saveRDS(pooled_results[[i]]$pooled_coef, file = paste0(path_for_pooled, "pooled_results_cluster", i, "_coef.rds"))
  saveRDS(pooled_results[[i]]$pooled_vcov, file = paste0(path_for_pooled, "pooled_results_cluster", i, "_vcov.rds"))
}

# Reconstruction
pred_reconstruct_pooled <- vector("list", num_clusters)
for (i in 1:num_clusters) {
  cluster_data <- read.csv("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/Non-imputed/final_temp_cluster_subgroup.csv")
  cluster_data <- as.data.table(cluster_data)
  cluster_data <- cluster_data %>% filter(!grepl("^2020", `year_month`))
  cluster_data <- cluster_data[cluster_data[['cluster_ward_std_6']] %in% c(i)] 
  
  Temp_measure <- "tmp_pw_percentile"
  minT <- min(cluster_data[[Temp_measure]], na.rm=TRUE)
  maxT <- max(cluster_data[[Temp_measure]], na.rm=TRUE)
  medT <- median(cluster_data[[Temp_measure]], na.rm=TRUE)
  quan01 <- quantile(cluster_data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
  quan25 <- quantile(cluster_data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
  quan75 <- quantile(cluster_data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
  quan99 <- quantile(cluster_data[[Temp_measure]], probs = 0.99, na.rm = TRUE)
  
  # Modeling (non-imputed; "original model")
  lagknots <- logknots(3, df = 3)
  knotstmean <- quantile(cluster_data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
  argvartmean <- list(fun="ns", knots=knotstmean)
  cbt <- crossbasis(cluster_data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=cluster_data$salid1)
  
  pred_reconstruct_pooled[[i]] <- crosspred(cbt, model.link = "log",
                                            coef = readRDS(paste0(path_for_pooled, "pooled_results_cluster", i, "_coef.rds")),
                                            vcov = readRDS(paste0(path_for_pooled, "pooled_results_cluster", i, "_vcov.rds")),
                                            cum = TRUE, cen = quan01)
}
