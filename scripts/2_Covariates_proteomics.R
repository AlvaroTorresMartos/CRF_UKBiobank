
# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, stringr, readr)

# Load the dataset
data <- readRDS("/mnt/project/data/processed/data_merged_copy_28_07_2025_needed_to_merge_by_ids_alvaro_2.rds")

# Remove metabolomics variables (columns starting with "metabo_")
metabo_vars <- grep("^metabo_", names(data), value = TRUE)
data <- data %>% select(-all_of(metabo_vars))

# Load proteomics data
proteomics <- read_csv("/mnt/project/data/raw/proteomics.csv")

# Merge datasets by "eid" (assuming 'eid' exists in both datasets)
data <- data %>% left_join(proteomics, by = "eid")

# Save the new dataset
saveRDS(data, "./upload/data_merged_copy_28_07_2025_needed_to_merge_by_ids_alvaro_2_proteomics.rds")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")