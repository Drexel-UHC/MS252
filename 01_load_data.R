library("tidyverse")
library("data.table")
library("dlnm"); library("gnm"); library("splines") # Modeling tools

# Load data
path_data <- "C:/Users/dpw48/OneDrive - Drexel University/git/MS252/Data/"

all <- fread(paste0(path_data, "final_temp_cluster_subgroup.csv"))