# Load required packages
if (!require("writexl")) install.packages("writexl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")
if (!require("purrr")) install.packages("purrr")
if (!require("readr")) install.packages("readr")
if (!require("gtsummary")) install.packages("gtsummary")

library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(writexl)
library(gtsummary)

# Load dataset
data <- readRDS("/mnt/project/data/processed/data_proteomics_10percent_filtered_25_07_2025.rds")


# Filter for non-missing VO2max
data_vo2 <- data %>% filter(is.na(`30038_0_0`))
n_total <- nrow(data_vo2)

# RDS
saveRDS(data_vo2, 
        "./upload/validation_cohort_proteomics.rds")

# Variable groups
continuous_vars <- c("30038_0_0")  # Add more continuous vars if needed
categorical_vars <- c("age_group", "bmi_group", "ethnicity", "townsend_q5",
                      "education_group", "region", "smoking_status", "pa_rec", "alc_group",
                      "fruit_veg_rec", "menopause", "sex", "sleep_group", "pregnancy",
                      "cancer", "cvd_fam", "diab_fam", "copd_fam")

# Variable labels
var_labels <- c(
  eid = "Participant ID",
  `30038_0_0` = "VO2max (ml/min/kg)",
  age_group = "Age group",
  bmi_group = "BMI group",
  ethnicity = "Ethnic background",
  townsend_q5 = "Townsend index (quintiles)",
  education_group = "Education level",
  region = "region",
  smoking_status = "Smoking status",
  pa_rec = "Meets PA recommendations",
  alc_group = "Alcohol intake frequency",
  fruit_veg_rec = "Meets fruit/veg recommendation",
  menopause = "Menopausal status",
  sex = "Sex",
  sleep_group = "Sleep duration group",
  pregnancy = "Pregnancy",
  cancer = "Cancer diagnosis",
  cvd_fam = "Family history: CVD",
  diab_fam = "Family history: Diabetes",
  copd_fam = "Family history: COPD"
)

# Table
summary_cont <- map_dfr(
  continuous_vars,
  function(var) {
    vec <- data_vo2[[var]]
    tibble(
      Variable = var_labels[[var]],
      Category = "All",
      N = sum(!is.na(vec)),
      `Percent of cohort` = round(100 * sum(!is.na(vec)) / n_total, 1),
      Mean = round(mean(vec, na.rm = TRUE), 2),
      SD = round(sd(vec, na.rm = TRUE), 2)
    )
  }
)

# Table
summary_cat <- map_dfr(
  categorical_vars,
  function(var) {
    df_summary <- data_vo2 %>%
      mutate(tmp_var = .data[[var]]) %>%
      group_by(tmp_var) %>%
      summarise(
        N = n(),
        Mean = mean(`30038_0_0`, na.rm = TRUE),
        SD = sd(`30038_0_0`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Variable = var_labels[[var]],
        Category = as.character(tmp_var)
      )
    
    df_summary <- df_summary %>%
      mutate(`Percent of cohort` = 100 * N / n_total)
    
    df_summary <- df_summary %>%
      mutate(`Percent of cohort` = round(`Percent of cohort`, 1))
    
    diff <- 100.0 - sum(df_summary$`Percent of cohort`)
    if (abs(diff) > 0) {
      idx <- which.max(df_summary$`Percent of cohort`)
      df_summary$`Percent of cohort`[idx] <- df_summary$`Percent of cohort`[idx] + diff
    }
    
    df_summary %>%
      mutate(
        Mean = round(Mean, 2),
        SD = round(SD, 2)
      ) %>%
      select(Variable, Category, N, `Percent of cohort`, Mean, SD)
  }
)

summary_final <- bind_rows(summary_cont, summary_cat)

write_xlsx(summary_final, "./upload/vo2max_summary_by_covariates_proteomics_validation_cohort.xlsx")

print(summary_final, n = 100)

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
