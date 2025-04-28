################################################################################
# Purpose: 
# Simulate demo datasets (5 cities x 60 days x 3 imputations), run DLNM for 
# main effects, subgroup analysis, and interaction effects. Apply Rubin's rule 
# and visualize pooled associations.
################################################################################

# ==== 1. Setup ====
library(data.table)
library(dlnm)
library(gnm)

set.seed(123)

n_cities <- 10
n_days <- 30
n_imputations <- 50

cities <- paste0("City", 1:n_cities)
dates <- seq.Date(from = as.Date("2020-01-01"), by = "day", length.out = n_days)

city_modifiers <- data.table(
  salid1 = cities,
  mean_temp_city = runif(n_cities, 20, 30),
  traffic_time = runif(n_cities, 30, 90),
  street_length = runif(n_cities, 5, 20)
)

# ==== 2. Simulate Imputed Datasets ====
simulate_imputation <- function(i) {
  expand.grid(salid1 = cities, allDate = dates) %>%
    as.data.table() %>%
    merge(city_modifiers, by = "salid1") %>%
    .[, `:=`(
      year = year(allDate),
      year_month_dow = paste0(format(allDate, "%Y-%m"), "-", weekdays(allDate)),
      tmp_pw_percentile = runif(.N, 0, 100),
      all_road = rpois(.N, lambda = 5),
      male_deaths = rbinom(.N, 5, 0.8),
      female_deaths = rbinom(.N, 5, 0.2)
    )]
}

demo_data_list <- lapply(1:n_imputations, simulate_imputation)

# ==== 3. Subgroup Analysis (Male Deaths Example) ====
coef_male <- list()
vcov_male <- list()

for (i in 1:n_imputations) {
  data <- demo_data_list[[i]]
  data[, stratum := factor(paste(salid1, year_month_dow))]
  
  cbt <- crossbasis(data$tmp_pw_percentile, lag = 2,
                    argvar = list(fun = "ns", knots = c(10, 75, 90)),
                    arglag = list(knots = logknots(3, df = 3)),
                    group = data$salid1)
  
  model_male <- gnm(male_deaths ~ cbt, eliminate = stratum, 
                    family = quasipoisson(), data = data, subset = male_deaths > 0)
  
  coef_male[[i]] <- coef(model_male)
  vcov_male[[i]] <- vcov(model_male)
}

# Rubin's Rule for Male Subgroup
coef_mat_male <- do.call(rbind, coef_male)
mean_coef_male <- colMeans(coef_mat_male)
var_within_male <- Reduce("+", vcov_male) / n_imputations
var_between_male <- cov(coef_mat_male)
total_vcov_male <- var_within_male + (1 + 1/n_imputations) * var_between_male

pred_male <- crosspred(cbt, model.link = "log", 
                       coef = mean_coef_male, vcov = total_vcov_male, 
                       cum = TRUE, cen = 50, by = 1)

# ==== 4. Interaction Analysis (Modifier: Mean City Temperature) ====
coef_interact1 <- list()
vcov_interact1 <- list()
coef_interact2 <- list()
vcov_interact2 <- list()

for (i in 1:n_imputations) {
  data <- demo_data_list[[i]]
  data[, stratum := factor(paste(salid1, year_month_dow))]
  
  # Crossbasis setup
  cbt <- crossbasis(data$tmp_pw_percentile, lag = 2,
                    argvar = list(fun = "ns", knots = c(25, 75)),
                    arglag = list(knots = logknots(3, df = 3)),
                    group = data$salid1)
  
  # Main model for interaction base
  model_road <- gnm(all_road ~ cbt, eliminate = stratum, 
                    family = quasipoisson(), data = data, subset = all_road > 0)
  
  # Define modifier and interaction terms
  modifier <- "mean_temp_city"
  intval <- quantile(data[[modifier]], c(0.10, 0.90))
  
  cbint1 <- cbt * (data[[modifier]] - intval[1])
  cbint2 <- cbt * (data[[modifier]] - intval[2])
  
  # Fit interaction models
  modint1 <- update(model_road, . ~ . + cbint1)
  modint2 <- update(model_road, . ~ . + cbint2)
  
  # Predictions
  pred_int1 <- crosspred(cbt, modint1, cen = 50)
  pred_int2 <- crosspred(cbt, modint2, cen = 50)
  
  # Store coefficients and vcov
  coef_interact1[[i]] <- coef(pred_int1)
  vcov_interact1[[i]] <- vcov(pred_int1)
  
  coef_interact2[[i]] <- coef(pred_int2)
  vcov_interact2[[i]] <- vcov(pred_int2)
}

# ==== Rubin's Rule for Interaction (10th Percentile) ====
coef_mat_int1 <- do.call(rbind, coef_interact1)
mean_coef_int1 <- colMeans(coef_mat_int1)
var_within_int1 <- Reduce("+", vcov_interact1) / n_imputations
var_between_int1 <- cov(coef_mat_int1)
total_vcov_int1 <- var_within_int1 + (1 + 1/n_imputations) * var_between_int1

pred_int1_pooled <- crosspred(cbt, model.link = "log",  
                              coef = mean_coef_int1, vcov = total_vcov_int1,  
                              cum = TRUE, cen = 50, by = 1)

# ==== Rubin's Rule for Interaction (90th Percentile) ====
coef_mat_int2 <- do.call(rbind, coef_interact2)
mean_coef_int2 <- colMeans(coef_mat_int2)
var_within_int2 <- Reduce("+", vcov_interact2) / n_imputations
var_between_int2 <- cov(coef_mat_int2)
total_vcov_int2 <- var_within_int2 + (1 + 1/n_imputations) * var_between_int2

pred_int2_pooled <- crosspred(cbt, model.link = "log",  
                              coef = mean_coef_int2, vcov = total_vcov_int2,  
                              cum = TRUE, cen = 50, by = 1)

# ==== 5. Plot Results ====
par(mfrow = c(1,3))
plot(pred_pooled, "overall", ylim = c(0.8, 1.6), main = "Main Effect", ylab = "RR")

plot(pred_male, "overall", ylim = c(0.8, 1.6), main = "Male Subgroup", ylab = "RR")

# Plot interaction comparison
plot(pred_int1_pooled, "overall", ylim = c(0.8, 1.6), 
     main = "Interaction: City Avg Temp", ylab = "RR", col = "blue")

lines(pred_int2_pooled, "overall", col = "red", lwd = 2)

legend("topleft", legend = c("10th percentile", "90th percentile"),
       col = c("blue", "red"), lty = 1, bty = "n")
