library("tidyverse")
library("data.table")
library("dlnm"); library("gnm"); library("splines") # Modeling tools

# Load data
path_data <- "C:/Users/dpw48/OneDrive - Drexel University/git/MS252/Data/"

all <- fread(paste0(path_data, "final_temp_cluster_subgroup.csv"))

# Select all temperature clusters for analysis
data <- all

# Select only specific temperature clusters for analysis
# data <- all[all[['cluster_ward_std_6']] %in% c(1)]