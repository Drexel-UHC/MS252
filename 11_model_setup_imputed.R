library("tidyverse")
library("data.table")
library("dlnm"); library("gnm"); library("splines") # Modeling tools


# Load data ---------------------------------------------------------------

# Number of imputations to run (max 100)
n_imp <- 5

# Load data
path_data <- "C:/Users/dpw48/git/MS252/Data/"

all <- fread(paste0(path_data, "final_temp_cluster_subgroup_1.csv"))

# Select temperature clusters for analysis (1-6, details in README)
temp_clusters <- 1:6

data <- all[cluster_ward_std_6 %in% temp_clusters]

# Set up model ------------------------------------------------------------

# Specify temperature measure
Temp_measure <- "tmp_pw_percentile" # percentile
# Temp_measure <- "L1ADtemp_pw"     # C°
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

# Specify knot number for lags
lagknots <- logknots(3, df = 3)

# Specify knot locations for temperature
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)

# Create crossbasis for DLNM
cbt <- crossbasis(data[[Temp_measure]], lag = 2,
                  argvar = argvartmean, arglag = list(knots = lagknots))

cbt_red <- crossbasis(data[[Temp_measure]], lag = 0,
                      argvar = argvartmean, group = data$salid1)


# Helper functions --------------------------------------------------------

# Accepts list of coefs and list cov matrices from analysis on imputed data
# Returns pooled estimates of the coefficients and variances via Rubin's rule
# https://documentation.sas.com/doc/en/statug/15.2/statug_mianalyze_details10.htm
apply_rubins_rule <- function(coef_list, vcov_list) {
  n_imp <- length(coef_list)
  
  # Get pooled estimate of coefs
  q_bar <- reduce(coef_list, `+`)/n_imp
  
  # Within-imputation variance
  W_bar <- reduce(vcov_list, `+`)/n_imp
  
  # Between-imputation covariance
  B <- cov(t(matrix(unlist(coef_list), ncol = n_imp)))
  
  # Total covariance
  T_0 <- W_bar + (1+1/n_imp)*B
  
  list(q_bar = q_bar,
       T_0 = T_0)
}