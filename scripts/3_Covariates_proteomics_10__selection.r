# COMPLEMENTARY ANALYSES OF PROTEOMIC DATA

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr)

# Load the dataset with integrated proteomics
data <- readRDS("/mnt/project/data/processed/data_merged_copy_28_07_2025_needed_to_merge_by_ids_alvaro_2_proteomics.rds")

# Load the original dataset without proteomics to identify which columns are new
original_data <- readRDS("/mnt/project/data/processed/data_merged_copy_28_07_2025_needed_to_merge_by_ids_alvaro_2.rds")

# Identify proteomics columns (those added after merging)
proteo_vars <- setdiff(names(data), names(original_data))

# Calculate proportion of missing values per subject across proteomics variables
na_proportion <- rowMeans(is.na(data[, proteo_vars]))

# Filter subjects with ≤10% missing values in proteomics
data_filtered <- data[na_proportion <= 0.10, ]

# Save the filtered dataset
saveRDS(data_filtered, "./upload/data_merged_copy_28_07_2025_needed_to_merge_by_ids_alvaro_2_proteomics_10_filtered.rds")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
