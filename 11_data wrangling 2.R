###########
# Set the directory containing the CSV files
data_directory <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS 252_all_imputed/Derived Data_processed"
output_directory <- "/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Data/MS 252_all_imputed/Derived Data_processed_by_imputation"

# Generate a dynamic list of road types from "road1" to "road100"
road_types <- paste0("road", 51:100)

# List all files in the directory
file_list <- list.files(path = data_directory, pattern = "\\.csv$", full.names = TRUE)

# Loop through each road type
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

###########