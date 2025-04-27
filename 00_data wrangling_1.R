# Load required libraries
library(haven)
library(lubridate)
library(RColorBrewer)
library(tsModel)
library(data.table)     # Efficient handling of large datasets
library(dlnm)           # Distributed lag nonlinear models
library(gnm)            # Generalized nonlinear models
library(splines)        # For spline functions
library(sf)             # Handle spatial data
library(terra)          # Handle raster data
library(exactextractr)  # Fast extraction of area-weighted raster cells
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stargazer)
library(broom)
library(tidyverse)
################################################################################
# Purpose:
# This script processes imputed road death data for 272 Latin American cities.
# Each city has a .sas7bdat file with 100 imputation columns (road1–road100),
# representing different sets of estimated daily road traffic deaths.
#
# The script performs the following:
# - Reads imputed SAS files for selected cities.
# - For each roadX imputation:
#     - Summarizes daily deaths overall and by subgroups (sex, age).
#     - Extracts and summarizes specific death types: vehicle, motorcycle,
#       bicycle, pedestrian.
#     - Merges each summary with daily temperature and city-level metadata.
#     - Adds calendar features (e.g., year, month, day of week).
#     - Saves each processed imputation summary to a separate CSV file.
#
# Input:  Imputed SAS datasets, city-level metadata (BEC, temperature clusters)
# Output: One processed CSV file per city and per roadX imputation
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################


#---------------------------------------------
# Step 1: Read in supporting city-level data
#---------------------------------------------
# Define directory and filename pattern
directory <- "/Volumes/TOSHIBA Kai/MS252/Derived Data_20240923"
file_pattern <- "^c(\\d+)\\.sas7bdat$"

# Read city-level modifier datasets
BEC1 <- fread("/Volumes/TOSHIBA Kai/MS252/MS252 requested data/BEC_L1AD_08162023.csv") %>%
  rename(salid1 = SALID1) %>%
  select(
    salid1, BECADSTTLGAVGL1AD, BECCZL1AD, BECADSTTDENSL1AD, BECADLRDENSL1AD, BECADINTDENSL1AD, BECPTCHDENSL1AD,
    BECADINTDENS3L1AD, BECADINTDENS4L1AD, BECADSTTPNODEAVGL1AD, BECADSTTPNODESDL1AD, BECADSTTLGAVGL1AD,
    BECADCRCTYAVGL1AD, BECSTTPL1AD, BECPCTURBANL1AD, BECGSPCTL1AD, BECGSPTCHDENSL1AD, BECMINWAGEL1AD,
    BECELEVATIONMAXL1AD, BECELEVATIOVEL1AD, BECELEVATIONMEDIANL1AD, BECELEVATIONMINL1AD, BECELEVATIONP25L1AD,
    BECELEVATIONP75L1AD, BECELEVATIONSTDL1AD, BECSLOPEMAXL1AD, BECSLOPEAVEL1AD, BECSLOPEMEDIANL1AD,
    BECSLOPEMINL1AD, BECSLOPEP25L1AD, BECSLOPEP75L1AD, BECSLOPESTDL1AD
  )

BEC2 <- fread("/Volumes/TOSHIBA Kai/MS252/MS252 requested data/BEC_RESTRICTED_L1AD_08162023.csv") %>%
  rename(salid1 = SALID1) %>%
  select(salid1, BECURBTRVDELAYINDEXL1AD, BECURBAVGTRAFTIMEL1AD, BECURBTRVDELAYTIMEL1AD, BECPARKPCTAREAL1AD)

TEMP_cluster <- read_sas("/Volumes/TOSHIBA Kai/MS252/city_level_temp_w_clusters.sas7bdat") %>%
  rename(salid1 = nsalid1) %>%
  select(salid1, cluster_ward_std_6, mean, std)

#---------------------------------------------
# Step 2: Load city-imputation files
#---------------------------------------------
# List all .sas7bdat files matching pattern
files <- list.files(path = directory, pattern = file_pattern, full.names = TRUE)

# Define subset for processing
half_point <- ceiling(length(files) / 2)
first_half_files <- files[1:half_point]
second_half_files <- files[(half_point + 1):length(files)]

# Example manual override to only process files from index 0 onwards
manual_files <- files[0:length(files)]

#---------------------------------------------
# Step 3: Define subgroups for summarization
#---------------------------------------------
subgroups <- list(
  list(name = "male", filter = "male == 1"),
  list(name = "female", filter = "male == 0"),
  list(name = "age1", filter = "age_cat == '<=9'"),
  list(name = "age2", filter = "age_cat == '10-19'"),
  list(name = "age3", filter = "age_cat == '20-34'"),
  list(name = "age4", filter = "age_cat == '35-64'"),
  list(name = "age5", filter = "age_cat == '65+'")
)

#---------------------------------------------
# Step 4: Define summarization function
#---------------------------------------------
summarize_by_group <- function(data, group_filter = NULL, death_var, result_name) {
  if (!is.null(group_filter)) {
    data <- data[eval(parse(text = group_filter))]
  }
  data %>%
    group_by(allDate) %>%
    summarise(!!result_name := sum(!!sym(death_var)), .groups = 'drop')
}

#---------------------------------------------
# Step 5: Process each city file
#---------------------------------------------
for (file_path in manual_files) {  # Replace with second_half_files if needed
  city_name <- sub("\\.sas7bdat$", "", basename(file_path))
  start_time <- Sys.time()
  cat("Processing city ID:", city_name, "\n")
  
  data <- as.data.table(read_sas(file_path))
  
  # Loop through each imputation (road1 to road100)
  for (road_num in 1:100) {
    road_var <- paste0("road", road_num)
    road_results <- list()
    
    # Total (unfiltered)
    result_name <- paste0("all_", road_var)
    road_results[[result_name]] <- summarize_by_group(data, NULL, road_var, result_name)
    
    # By subgroups
    for (subgroup in subgroups) {
      result_name <- paste0(subgroup$name, "_", road_var)
      road_results[[result_name]] <- summarize_by_group(data, subgroup$filter, road_var, result_name)
    }
    
    # Additional death types
    other_death_types <- c("vehicle", "motorcycle", "bicycle", "ped")
    for (death_type in other_death_types) {
      death_var <- paste0(death_type, road_num)
      result_name <- paste0(death_type, "_", road_var)
      road_results[[result_name]] <- summarize_by_group(data, NULL, death_var, result_name)
    }
    
    # Combine all summaries
    road_final_result <- Reduce(function(x, y) full_join(x, y, by = "allDate"), road_results)
    
    # Add temporal variables
    road_final_result <- road_final_result %>%
      mutate(
        allDate = as.Date(allDate),
        year = year(allDate),
        month = lubridate::month(allDate, label = TRUE),
        dow = lubridate::wday(allDate, label = TRUE),
        year_month = paste0(year, "-", month),
        year_month_dow = paste0(year, "-", month, "-", dow),
        day_of_year = yday(allDate)
      )
    
    # Extract city/date-level attributes
    constant_columns <- data %>%
      select(allDate, country, salid1) %>%
      distinct()
    
    varying_columns <- data %>%
      select(allDate, L1ADtemp_pw, L1ADtemp_x, tmp_pw_percentile, tmp_x_percentile, city_size) %>%
      distinct()
    
    # Merge all data
    road_final_result <- road_final_result %>%
      left_join(constant_columns, by = "allDate") %>%
      left_join(varying_columns, by = "allDate") %>%
      left_join(BEC1, by = "salid1") %>%
      left_join(BEC2, by = "salid1") %>%
      left_join(TEMP_cluster, by = "salid1")
    
    # Export CSV
    out_directory <- "/Volumes/TOSHIBA Kai/MS252/Derived Data_20240923_processed"
    output_file_path <- file.path(out_directory, paste0(city_name, "_road", road_num, "_processed.csv"))
    fwrite(road_final_result, output_file_path)
    
    # Clean up memory
    rm(road_results, road_final_result)
    gc()
  }
  
  end_time <- Sys.time()
  cat("Time taken for city ID:", city_name, "is", end_time - start_time, "\n")
  
  rm(data, constant_columns, varying_columns)
  gc()
}
