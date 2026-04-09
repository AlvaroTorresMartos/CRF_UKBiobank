# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026
# Table 1 - Metabolomic validation cohort (excluding discovery cohort)
# NOTE: This cohort has NO VO2max -> frequency table only

# ---------------------------
# 1) Packages
# ---------------------------
pkg_needed <- c(
  "dplyr", "tibble", "purrr", "readr",
  "officer", "flextable", "stringr"
)

to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(officer)
library(flextable)
library(stringr)

# ---------------------------
# 2) Load dataset (base completa)
# ---------------------------
data <- readRDS(
  "/mnt/project/data/processed/metabolomics_clean_without_NCDs_17022026.rds"
)

# ---------------------------
# 3) Filter for cohort with metabolomics 
#    BUT excluding individuals from discovery cohort
# ---------------------------
discovery_cohort <- readRDS(
  "/mnt/project/data/processed/discovery_cohort_metabolomics.rds"
)

eid_discovery <- discovery_cohort$eid

data_validation <- data %>%
  filter(!eid %in% eid_discovery)

cat("Validation cohort size (raw):", nrow(data_validation), "\n")

# Save validation cohort
saveRDS(
  data_validation,
  "./upload/validation_cohort_metabolomics.rds"
)

# ---------------------------
# 4) Variable groups
# ---------------------------
categorical_vars <- c(
  "age_group", "bmi_group", "ethnicity", "townsend_q5",
  "education_group", "region", "smoking_status", "pa_rec", "alc_group",
  "fruit_veg_rec", "menopause", "sex", "sleep_group", "pregnancy",
  "cancer", "cvd_fam", "diab_fam", "copd_fam"
)

# ---------------------------
# 5) Variable labels (match manuscript wording)
# ---------------------------
var_labels <- c(
  eid = "Participant ID",
  
  age_group = "Age (years)",
  bmi_group = "Body mass index (kg·m-2)",
  ethnicity = "Ethnic background",
  townsend_q5 = "Townsend deprivation index",
  education_group = "Qualifications",
  region = "Geographical region",
  smoking_status = "Smoking status",
  pa_rec = "Meet physical activity recommendations",
  alc_group = "Alcohol intake",
  fruit_veg_rec = "Meet recommended fruit and vegetable intake",
  menopause = "History of menopause",
  sex = "Sex",
  sleep_group = "Sleep duration group",
  pregnancy = "Pregnancy",
  cancer = "Family history of cancer",
  cvd_fam = "Family history of cardiovascular disease",
  diab_fam = "Family history of diabetes mellitus",
  copd_fam = "Family history: COPD"
)

# ---------------------------
# 6) Sanity checks + build dataset for table
# ---------------------------
missing_cat <- setdiff(categorical_vars, names(data_validation))
if (length(missing_cat) > 0) {
  warning("⚠️ Missing categorical variables in data_validation: ",
          paste(missing_cat, collapse = ", "),
          "\nThey will be skipped.")
  categorical_vars <- intersect(categorical_vars, names(data_validation))
}

# Build dataset: only eid + categorical variables
# Convert NAs to "Unknown" so they appear explicitly in Table 1 (optional but recommended)
data_tab <- data_validation %>%
  select(eid, all_of(categorical_vars)) %>%
  mutate(across(all_of(categorical_vars), ~ ifelse(is.na(.x) | .x == "", "Unknown", as.character(.x))))

n_total <- nrow(data_tab)
cat("Validation cohort size (used for table):", n_total, "\n")

# ---------------------------
# 7) Formatting helpers
# ---------------------------
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct1 <- function(x) sprintf("%.1f", x)

# ---------------------------
# 8) Build table blocks (paper style)
#     - Characteristic shown only on first row of each block
# ---------------------------
build_block <- function(var, label, n_total, df) {
  
  out <- df %>%
    mutate(.grp = .data[[var]]) %>%
    group_by(.grp) %>%
    summarise(
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(
      pct  = 100 * N / n_total,
      freq = sprintf("%s (%s%%)", fmt_int(N), fmt_pct1(pct)),
      Characteristic = label,
      Category = as.character(.grp)
    ) %>%
    select(
      Characteristic,
      Category,
      `Frequency of participants (% full cohort)` = freq
    )
  
  # Show characteristic label only once per block (Table 1 style)
  if (nrow(out) > 1) out$Characteristic[-1] <- ""
  
  out
}

# ---------------------------
# 9) Create final table (top "All" row + blocks)
# ---------------------------
all_row <- tibble(
  Characteristic = "",
  Category = "",
  `Frequency of participants (% full cohort)` = sprintf("%s (100.0%%)", fmt_int(n_total))
)

cat_blocks <- map_dfr(
  categorical_vars,
  ~ build_block(.x, var_labels[[.x]], n_total, data_tab)
)

final_table <- bind_rows(all_row, cat_blocks)

# ---------------------------
# 10) Make flextable and Word export
# ---------------------------
ft <- flextable(final_table)

ft <- set_header_labels(
  ft,
  Characteristic = "Characteristic",
  Category = "Category",
  `Frequency of participants (% full cohort)` = "Frequency of participants (% full cohort)"
)

# Styling: clean, “paper-ish”
ft <- theme_vanilla(ft)
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 10, part = "all")
ft <- align(ft, align = "left", part = "all")
ft <- align(
  ft,
  j = "Frequency of participants (% full cohort)",
  align = "center",
  part = "all"
)
ft <- padding(ft, padding = 2, part = "all")
ft <- autofit(ft)

# Borders: subtle journal-like rules
ft <- border_remove(ft)
std_border <- fp_border(color = "black", width = 0.75)
ft <- hline_top(ft, border = std_border, part = "all")
ft <- hline(ft, border = fp_border(color = "black", width = 0.25), part = "body")
ft <- hline_bottom(ft, border = std_border, part = "all")

# ---------------------------
# 11) Export to Word
# ---------------------------
out_docx <- "./upload/validation_cohort_metabolomics.docx"

doc <- read_docx()
doc <- body_add_par(doc, "METABOLOMIC VALIDATION COHORT", style = "heading 1")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_flextable(doc, value = ft)

print(doc, target = out_docx)

message("✅ Saved Word table to: ", out_docx)

# ---------------------------
# 12) Print preview in console (optional)
# ---------------------------
print(final_table, n = 200)

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
