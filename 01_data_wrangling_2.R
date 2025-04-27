################################################################################
# Purpose:
# This script aggregates previously processed city-imputation CSV files into
# imputation-level datasets. Specifically:
#
# Input:
#   - 27,200 CSV files (272 cities × 100 imputations), each named like:
#       [city_id]_road[1–100]_processed.csv
#
# Output:
#   - 100 CSV files, each representing a single imputation across all 272 cities,
#     named: road1.csv, road2.csv, ..., road100.csv
#
# Example use case:
#   - Enables imputation-level modeling by combining all cities for a given
#     imputation round.
#
# Author: Cheng-Kai Hsu
# Date: April 2025
################################################################################

#---------------------------------------------
# Step 1: Define directories and road types
#---------------------------------------------
data_directory <- "/Volumes/TOSHIBA Kai/MS252/Derived Data_20240923_processed"
output_directory <- "/Volumes/TOSHIBA Kai/MS252/Derived Data_20240923_processed_by_imputation"

# Define which imputations to process (adjust range as needed)
# Generate a dynamic list of road types e.g., from "road1" to "road100"
road_types <- paste0("road", 21:100)

#---------------------------------------------
# Step 2: List all available CSV files
#---------------------------------------------
file_list <- list.files(path = data_directory, pattern = "\\.csv$", full.names = TRUE)

#---------------------------------------------
# Step 3: Aggregate files by imputation type
#---------------------------------------------
for (road_type in road_types) {
  # Record the start time
  start_time <- Sys.time()  
  
  # Create a pattern to match the current road type in file names, ensuring an exact match
  pattern <- paste0("_", road_type, "_")
  
  # Filter files based on the generated pattern for the current road type
  filtered_files <- file_list[grepl(pattern, file_list)]
  
  # Read all filtered files for the current road type into a list of data frames
  road_data_list <- lapply(filtered_files, read.csv)
  
  # Combine all data frames into a single data frame for the current road type
  combined_data <- do.call(rbind, road_data_list) 
  
  # Create the output file name
  output_file <- file.path(output_directory, paste0(road_type, ".csv"))
  
  # Write the combined data to a CSV file
  write.csv(combined_data, file = output_file, row.names = FALSE)
  
  # Record the end time
  end_time <- Sys.time()
  # Calculate and print the elapsed time
  elapsed_time <- end_time - start_time
  cat("Time taken for", road_type, "is", elapsed_time, "\n")  
}
