################################################################################
# Purpose:
# This script processes *median-imputed* road death data for 272 Latin American cities.
# Each city originally contains 100 imputed estimates of daily road traffic deaths (road1–road100),
# but this script uses only the *median* across those 100 imputations as the outcome variable.
#
# The script performs the following:
# - Loads summarized datasets where the median of 100 imputations has been precomputed.
# - Aggregates daily deaths overall and by subgroups (sex, age).
# - Adds calendar features (e.g., year, month, day of week).
# - Merges summaries with temperature data and city-level metadata (BEC, SEC, REG).
# - Fits distributed lag non-linear models (DLNMs) using temperature as the exposure.
#
# Input:  Median-imputed road death datasets, city-level metadata (BEC, temperature clusters)
# Output: Cleaned summary datasets and fitted DLNM models
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

# Load required libraries
library(tidyverse)
library(haven)
library(lubridate)
library(data.table)
library(dlnm)
library(gnm)

# 1. Paths to Data ----------------------------------------------------------
path_1 <- "/Volumes/TOSHIBA Kai/MS252/Median-imputed/Prelim Derived Data_20240920/"
path_write <- "/Volumes/TOSHIBA Kai/MS252/Median-imputed/"
path_2 <- "/Volumes/TOSHIBA Kai/MS252/MS252 requested data/"
path_3 <- "/Volumes/TOSHIBA Kai/MS252"

# 2. Helper Functions -------------------------------------------------------

# Function to add date-related variables to dataset
add_year_month_dow <- function(df) {
  df %>% 
    mutate(
      age_cat = factor(age_cat),
      year = str_extract(allDate, "\\d{4}"),
      year_month = str_extract(allDate, "\\d{4}-\\d{2}"),
      doy = yday(allDate),
      dow = lubridate::wday(allDate, label = TRUE),
      year_month_dow = paste(year_month, dow, sep = "-")
    )
}

# Custom summarization function to aggregate deaths and median values
custom_summary_function <- function(df) {
  df %>% 
    summarize(
      across(c(country, city_size, salid1, year_month, dow), first),
      across(c(deaths, median_road_round, median_vehicle_round, median_motorcycle_round,
               median_bicycle_round, median_ped_round), sum),
      across(c(L1ADtemp_pw, L1ADtemp_x, tmp_pw_percentile, tmp_x_percentile), mean),
      male_deaths    = sum(median_road * (male == 1)),
      female_deaths  = sum(median_road * (male == 0)),
      deaths_under_9 = sum(median_road * (age_cat == "<=9")),
      deaths_10_19   = sum(median_road * (age_cat == "10-19")),
      deaths_20_34   = sum(median_road * (age_cat == "20-34")),
      deaths_35_64   = sum(median_road * (age_cat == "35-64")),
      deaths_65_plus = sum(median_road * (age_cat == "65+"))
    )
}

# 3. Load and Summarize Mortality Data ---------------------------------------
data_files_1 <- list.files(path_1, pattern = ".sas7bdat")

# Read and summarize data
combined_df <- data_files_1 %>% 
  map_dfr(\(x) {
    haven::read_sas(file.path(path_1, x)) %>% 
      add_year_month_dow() %>% 
      group_by(allDate) %>% 
      custom_summary_function()
  })

# Save combined dataset
# write_csv(combined_df, file.path(path_write, "combined_subgroup.csv"))

# 4. Load City-Level Modifier Data -------------------------------------------

# Load BEC data
BEC1 <- read_csv(file.path(path_2, "BEC_L1AD_08162023.csv")) %>%
  select(matches("ISO2|SALID1|BEC.*"))

BEC2 <- read_csv(file.path(path_2, "BEC_RESTRICTED_L1AD_08162023.csv")) %>%
  select(matches("SALID1|BEC.*"))

BEC <- inner_join(BEC1, BEC2, by = "SALID1")

# Load SEC data
SEC <- read_csv(file.path(path_2, 'SEC_INDEXSCORES_L1AD_07102023.csv')) %>% 
  filter(YEAR == max(YEAR), .by = SALID1) %>%
  select(SALID1, CNSSEI_L1AD, CNSSE1_L1AD, CNSSE2_L1AD, CNSSE3_L1AD)

# Load REG data
REG <- read_csv(file.path(path_2, 'VehicleRegistration_L1AD_20201027.csv')) %>% 
  filter(YEAR == max(YEAR), .by = SALID1) %>%
  select(SALID1, BECMTRBRATEL1AD, BECPAVRATEL1AD, BECTOTVRRATEL1AD)

# Combine mortality and modifiers
df_all <- combined_df %>% 
  left_join(BEC, by = c("salid1" = "SALID1")) %>% 
  left_join(SEC, by = c("salid1" = "SALID1")) %>% 
  left_join(REG, by = c("salid1" = "SALID1"))

# Save combined data
# write_csv(df_all, file.path(path_write, "final_subgroup.csv"))

# 5. Merge Temperature Cluster Data -----------------------------------------
temp_cluster <- haven::read_sas(file.path(path_3, "city_level_temp_w_clusters.sas7bdat"))

df_temp_cluster <- df_all %>% 
  left_join(temp_cluster, by = c("salid1" = "nsalid1"))

# Save final merged data
# write_csv(df_temp_cluster, file.path(path_write, "final_temp_cluster_subgroup.csv"))

# 6. Load Final Dataset -----------------------------------------------------
median_imputed <- fread(file.path(path_write, "final_temp_cluster_subgroup.csv"))

median_imputed <- median_imputed %>%
  mutate(year = lubridate::year(as.Date(allDate))) %>%
  filter(!(salid1 == 204106 & year >= 2011 & year <= 2020))

# Exclude 2020 and prepare strata
data <- median_imputed %>% filter(!grepl("^2020", year_month)) %>% as.data.table()
rm(median_imputed)
data[, stratum := factor(paste(salid1, year_month, dow, sep = ":"))]
data[, keep := sum(median_road_round) > 0, by = stratum]

# 7. Temperature Variables Setup --------------------------------------------
Temp_measure <- "tmp_pw_percentile"
Temp_measure_celcius <- "L1ADtemp_pw"

# Calculate quantiles (percentiles)
minT <- min(data[[Temp_measure]], na.rm=TRUE)
maxT <- max(data[[Temp_measure]], na.rm=TRUE)
medT <- median(data[[Temp_measure]], na.rm=TRUE)
quan01 <- quantile(data[[Temp_measure]], probs = 0.01, na.rm = TRUE)
quan25 <- quantile(data[[Temp_measure]], probs = 0.25, na.rm = TRUE)
quan75 <- quantile(data[[Temp_measure]], probs = 0.75, na.rm = TRUE)
quan99 <- quantile(data[[Temp_measure]], probs = 0.99, na.rm = TRUE)

Temp_measure_celcius <- "L1ADtemp_pw"
minT_celcius <- min(data[[Temp_measure_celcius]], na.rm=TRUE)
maxT_celcius <- max(data[[Temp_measure_celcius]], na.rm=TRUE)
medT_celcius <- median(data[[Temp_measure_celcius]], na.rm=TRUE)
quan01_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.01, na.rm = TRUE)
quan25_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.25, na.rm = TRUE)
quan75_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.75, na.rm = TRUE)
quan99_celcius <- quantile(data[[Temp_measure_celcius]], probs = 0.99, na.rm = TRUE)


# 8. Modeling ---------------------------------------------------------------
lagknots <- logknots(3, df = 3)
knotstmean <- quantile(data[[Temp_measure]], c(10,75,90)/100, na.rm=T)
argvartmean <- list(fun="ns", knots=knotstmean)
cbt <- crossbasis(data[[Temp_measure]], lag = 2, argvar = argvartmean, arglag = list(knots = lagknots), group=data$salid1)
model_road <- gnm(median_road_round ~ cbt, eliminate = stratum, 
                  family = quasipoisson(), data = data, subset=keep)

lagknots_celcius <- logknots(3, df = 3)
knotstmean_celcius <- quantile(data[[Temp_measure_celcius]], c(10,75,90)/100, na.rm=T)
argvartmean_celcius <- list(fun="ns", knots=knotstmean_celcius)
cbt_celcius <- crossbasis(data[[Temp_measure_celcius]], lag = 2, argvar = argvartmean_celcius, arglag = list(knots = lagknots_celcius), group=data$salid1)
model_road_celcius <- gnm(median_road_round ~ cbt_celcius, eliminate = stratum, 
                          family = quasipoisson(), data = data, subset=keep)

pred_road <- crosspred(cbt, model_road, cum=TRUE, cen=quan01, by=0.1)
pred_road_celcius <- crosspred(cbt_celcius, model_road_celcius, cum=TRUE, cen=quan01_celcius)

