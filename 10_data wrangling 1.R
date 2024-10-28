library(haven)
library(lubridate)
library(RColorBrewer)
library(tsModel)
library(data.table) # HANDLE LARGE DATASETS
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(sf) ; library(terra) # HANDLE SPATIAL DATA
library(exactextractr) # FAST EXTRACTION OF AREA-WEIGHTED RASTER CELLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
library(stargazer)
library(broom)
library("tidyverse")
library("lubridate")

######This script turns the raw data (272 files, each representing a city with 100 imputations) in to 27200 files (each representing a city-imputation)

# Define the directory and file pattern
directory <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Derived Data_20240923/"
file_pattern <- "^c(\\d+)\\.sas7bdat$"

# Read city-level modifiers
BEC1 <- fread('/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252 requested data/BEC_L1AD_08162023.csv')
BEC2 <- fread('/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252 requested data/BEC_RESTRICTED_L1AD_08162023.csv')
TEMP_cluster <- read_sas('/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/city_level_temp_w_clusters.sas7bdat')

BEC1 <- BEC1 %>%
  rename(salid1 = SALID1) %>%
  select(salid1, BECADSTTLGAVGL1AD,BECCZL1AD,BECADSTTDENSL1AD,BECADLRDENSL1AD,BECADINTDENSL1AD,BECPTCHDENSL1AD,BECADINTDENS3L1AD,BECADINTDENS4L1AD,BECADSTTPNODEAVGL1AD,BECADSTTPNODESDL1AD,BECADSTTLGAVGL1AD,BECADCRCTYAVGL1AD,BECSTTPL1AD,BECPCTURBANL1AD,BECGSPCTL1AD,BECGSPTCHDENSL1AD,BECMINWAGEL1AD,BECELEVATIONMAXL1AD,BECELEVATIOVEL1AD,BECELEVATIONMEDIANL1AD,BECELEVATIONMINL1AD,BECELEVATIONP25L1AD,BECELEVATIONP75L1AD,BECELEVATIONSTDL1AD,BECSLOPEMAXL1AD,BECSLOPEAVEL1AD,BECSLOPEMEDIANL1AD,BECSLOPEMINL1AD,BECSLOPEP25L1AD,BECSLOPEP75L1AD,BECSLOPESTDL1AD)

BEC2 <- BEC2 %>%
  rename(salid1 = SALID1) %>%
  select(salid1, BECURBTRVDELAYINDEXL1AD,BECURBAVGTRAFTIMEL1AD,BECURBTRVDELAYTIMEL1AD,BECPARKPCTAREAL1AD)

TEMP_cluster <- TEMP_cluster %>%
  rename(salid1 = nsalid1) %>%
  select(salid1, cluster_ward_std_6, mean, std) 

# List all files in the specified directory with the specified pattern
files <- list.files(path = directory, pattern = file_pattern, full.names = TRUE)

# Split files into two halves
half_point <- ceiling(length(files) / 2)
first_half_files <- files[1:half_point]  # Process first half
second_half_files <- files[(half_point + 1):length(files)]  # Uncomment for the second half

manual_files <- files[(255):length(files)]

# Define the subgroups and death types to summarize
subgroups <- list(
  list(name = "male", filter = "male == 1"),
  list(name = "female", filter = "male == 0"),
  list(name = "age1", filter = "age_cat == '<=9'"),
  list(name = "age2", filter = "age_cat == '10-19'"),
  list(name = "age3", filter = "age_cat == '20-34'"),
  list(name = "age4", filter = "age_cat == '35-64'"),
  list(name = "age5", filter = "age_cat == '65+'")
)

# Function to summarize road deaths by subgroup for a specific death variable (e.g., road1, road2)
summarize_by_group <- function(data, group_filter = NULL, death_var, result_name) {
  if (!is.null(group_filter)) {
    data <- data[eval(parse(text = group_filter))]
  }
  
  data %>%
    group_by(allDate) %>%
    summarise(!!result_name := sum(!!sym(death_var)), .groups = 'drop')
}

# Loop through each file (city) in the first half
for (file_path in manual_files) { #replace "first_half_files" to "second_half_files" for second batch
  # Extract city identifier from the file name
  city_name <- sub("\\.sas7bdat$", "", basename(file_path))
  
  # Record the start time
  start_time <- Sys.time()
  
  # Print the current city being processed
  cat("Processing city ID:", city_name, "\n")
  
  # Read the SAS file for the current city
  data <- as.data.table(read_sas(file_path))
  
  # Loop over road types to create summaries for road and other death types
  for (road_num in 1:100) {
    road_var <- paste0("road", road_num)
    
    # Initialize an empty list to store results for the current road type
    road_results <- list()
    
    # Add the all-road summary without filtering
    result_name <- paste0("all_", road_var)
    road_results[[result_name]] <- summarize_by_group(data, NULL, road_var, result_name)
    
    # Summarize for each subgroup
    for (subgroup in subgroups) {
      result_name <- paste0(subgroup$name, "_", road_var)
      road_results[[result_name]] <- summarize_by_group(data, subgroup$filter, road_var, result_name)
    }
    
    # Summarize other death types (vehicle, motorcycle, bicycle, pedestrian) for each road type
    other_death_types <- c("vehicle", "motorcycle", "bicycle", "ped")
    for (death_type in other_death_types) {
      death_var <- paste0(death_type, road_num) # Create death variable name dynamically
      result_name <- paste0(death_type, "_", road_var)
      road_results[[result_name]] <- summarize_by_group(data, NULL, death_var, result_name)
    }
    
    # Combine all results for the current road type into one data frame
    road_final_result <- Reduce(function(x, y) full_join(x, y, by = "allDate"), road_results)
    
    # Ensure 'allDate' is in Date format and create additional time-related columns
    road_final_result <- road_final_result %>%
      mutate(
        allDate = as.Date(allDate),
        year = year(allDate),
        month = month(allDate, label = TRUE),
        dow = wday(allDate, label = TRUE),
        year_month = paste0(year, "-", month),
        year_month_dow = paste0(year, "-", month, "-", dow),
        day_of_year = yday(allDate)
      )
    
    # Extract constant and varying columns from the original data
    constant_columns <- data %>%
      select(allDate, country, salid1) %>%
      distinct()
    
    varying_columns <- data %>%
      select(allDate, L1ADtemp_pw, L1ADtemp_x, tmp_pw_percentile, tmp_x_percentile, city_size) %>%
      distinct()
    
    # Merge constant and varying columns with the road-specific final result
    road_final_result <- road_final_result %>%
      left_join(constant_columns, by = "allDate") %>%
      left_join(varying_columns, by = "allDate")
    
    # Merge additional datasets on 'salid1'
    road_final_result <- road_final_result %>%
      left_join(BEC1, by = "salid1") %>%
      left_join(BEC2, by = "salid1") %>%
      left_join(TEMP_cluster, by = "salid1")
    
    # Save the road-specific result to a CSV file
    out_directory <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS252_impandnonimp_Sep24/imputed/Derived Data_20240923_processed"
    output_file_path <- file.path(out_directory, paste0(city_name, "_road", road_num, "_processed.csv"))
    fwrite(road_final_result, output_file_path)
    
    # Clean up memory after processing each road type
    rm(road_results, road_final_result)
    gc()  # Garbage collection to free up memory
  }
  
  # Record the end time
  end_time <- Sys.time()
  # Calculate and print the elapsed time
  elapsed_time <- end_time - start_time
  cat("Time taken for city ID:", city_name, "is", elapsed_time, "\n")
  
  # Clean up memory after processing each city
  rm(data, constant_columns, varying_columns)
  gc()  # Garbage collection to free up memory
}
