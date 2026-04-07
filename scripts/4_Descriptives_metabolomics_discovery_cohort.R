# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026

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
# 2) Load dataset
# ---------------------------
data <- readRDS(
  "/metabolomics_clean_without_NCDs_17022026.rds"
)

# ---------------------------
# 3) Filter for non-missing VO2max
# ---------------------------
vo2_var <- "30038_0_0"  # VO2max field
data_vo2 <- data %>% filter(!is.na(.data[[vo2_var]]))
n_total <- nrow(data_vo2)

# Save filtered subcohort (RDS)
saveRDS(
  data_vo2,
  "/Final databases/discovery_cohort_metabolomics.rds"
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
# 5) Variable labels (edit to match manuscript wording exactly)
# ---------------------------
var_labels <- c(
  eid = "Participant ID",
  `30038_0_0` = "VO2max (ml/min/kg)",
  
  # If you want EXACT paper labels (like your example), I recommend:
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
# 6) Formatting helpers
# ---------------------------
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct1 <- function(x) sprintf("%.1f", x)
fmt_mean_sd <- function(m, s, digits = 1) sprintf(paste0("%.", digits, "f \u00B1 %.", digits, "f"), m, s)

# ---------------------------
# 7) Build table blocks (paper style)
#     - Characteristic shown only on first row of each block
# ---------------------------
build_block <- function(var, label, vo2_var, n_total, df) {
  
  # Make sure missing categories aren't silently dropped (optional):
  # df <- df %>% mutate("{var}" := ifelse(is.na(.data[[var]]), "Unknown", as.character(.data[[var]])))
  
  out <- df %>%
    mutate(.grp = .data[[var]]) %>%
    group_by(.grp) %>%
    summarise(
      N = n(),
      Mean = mean(.data[[vo2_var]], na.rm = TRUE),
      SD   = sd(.data[[vo2_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pct  = 100 * N / n_total,
      freq = sprintf("%s (%s%%)", fmt_int(N), fmt_pct1(pct)),
      crf  = fmt_mean_sd(Mean, SD, digits = 1),
      Characteristic = label,
      Category = as.character(.grp)
    ) %>%
    select(
      Characteristic,
      Category,
      `Frequency of participants (% full cohort)` = freq,
      `CRF (VO2max; ml/min/kg)` = crf
    )
  
  # Show characteristic label only once per block (classic Table 1 style)
  if (nrow(out) > 1) out$Characteristic[-1] <- ""
  
  out
}

# ---------------------------
# 8) Create final table (top "All" row + blocks)
# ---------------------------
all_row <- tibble(
  Characteristic = "",
  Category = "",
  `Frequency of participants (% full cohort)` = sprintf("%s (100.0%%)", fmt_int(n_total)),
  `CRF (VO2max; ml/min/kg)` = fmt_mean_sd(
    mean(data_vo2[[vo2_var]], na.rm = TRUE),
    sd(data_vo2[[vo2_var]], na.rm = TRUE),
    digits = 1
  )
)

cat_blocks <- map_dfr(
  categorical_vars,
  ~ build_block(.x, var_labels[[.x]], vo2_var, n_total, data_vo2)
)

final_table <- bind_rows(all_row, cat_blocks)

# ---------------------------
# 9) Make flextable with requested schema and Word export
# ---------------------------
ft <- flextable(final_table)

ft <- set_header_labels(
  ft,
  Characteristic = "Characteristic",
  Category = "Category",
  `Frequency of participants (% full cohort)` = "Frequency of participants (% full cohort)",
  `CRF (VO2max; ml/min/kg)` = "CRF\n(VO2max; ml/min/kg)"
)

# Styling: clean, “paper-ish”
ft <- theme_vanilla(ft)
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 10, part = "all")
ft <- align(ft, align = "left", part = "all")
ft <- align(
  ft,
  j = c("Frequency of participants (% full cohort)", "CRF (VO2max; ml/min/kg)"),
  align = "center",
  part = "all"
)
ft <- padding(ft, padding = 2, part = "all")
ft <- autofit(ft)

# Optional: reduce visual “Excel vibe” (more journal-like)
# Remove most borders and keep subtle horizontal rules
ft <- border_remove(ft)
std_border <- fp_border(color = "black", width = 0.75)
ft <- hline_top(ft, border = std_border, part = "all")
ft <- hline(ft, border = fp_border(color = "black", width = 0.25), part = "body")
ft <- hline_bottom(ft, border = std_border, part = "all")

# ---------------------------
# 10) Export to Word
# ---------------------------
out_docx <- "/discovery_cohort_metabolomics.docx"

doc <- read_docx()
doc <- body_add_par(doc, "METABOLOMIC DISCOVERY COHORT", style = "heading 1")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_flextable(doc, value = ft)

print(doc, target = out_docx)

message("✅ Saved Word table to: ", out_docx)

# ---------------------------
# 11) Print preview in console (optional)
# ---------------------------
print(final_table, n = 200)
