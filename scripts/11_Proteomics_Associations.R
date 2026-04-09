
# ==========================
# 1. LOAD REQUIRED PACKAGES
# ==========================
packages_needed <- c("broom", "dplyr", "lm.beta", "writexl")

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(broom)
library(dplyr)
library(lm.beta)
library(writexl)

# ==========================
# 2. MAIN FUNCTION
# ==========================
recursive_lm <- function(outcome, predictor, data, formula, min_n = 50) {
  
  results <- lapply(outcome, function(y) {
    
    covariates <- trimws(unlist(strsplit(formula, "\\+")))
    covariates <- covariates[covariates != ""]
    
    model_vars <- c(y, predictor, covariates)
    data_model <- data[complete.cases(data[, model_vars]), ]
    
    cat("Protein:", y, "- N =", nrow(data_model), "\n")
    
    if (nrow(data_model) < min_n) return(NULL)
    
    object2 <- lm(
      as.formula(paste0(y, " ~ `", predictor, "`", formula)), data = data_model
    )
    
    tidy_object2 <- broom::tidy(object2, conf.int = TRUE) %>%
      filter(term == predictor | term == paste0("`", predictor, "`"))
    
    if (nrow(tidy_object2) == 0) {
      cat("Predictor term not found in:", y, "\n")
      return(NULL)
    }
    
    object1 <- lm.beta::lm.beta(object2)
    
    tidy_object1 <- broom::tidy(object1) %>%
      select(term, std_estimate) %>%
      rename(predictor = term, std_coefficient = std_estimate)
    
    tidy_object2 <- tidy_object2 %>%
      rename(predictor = term) %>%
      mutate(
        outcome = y,
        formula = paste(attr(object2$terms, "term.labels"), collapse = " "),
        N = nrow(data_model)
      ) %>%
      relocate(outcome, .before = predictor) %>%
      rename(
        beta_coefficient = estimate,
        std_error = std.error,
        t_statistic = statistic,
        CI_low = conf.low,
        CI_high = conf.high
      )
    
    tidy_object <- merge(tidy_object1, tidy_object2, by = "predictor") %>%
      relocate(outcome, .before = predictor) %>%
      relocate(beta_coefficient, .before = std_coefficient)
    
    return(tidy_object)
  })
  
  final_result <- bind_rows(results)
  return(final_result)
}

# ==========================
# 3. LOAD DATA & DEFINE VARIABLES (PROTEINS)
# ==========================
data_merged_filtered_proteomics <- readRDS(
  "/mnt/project/data/processed/data_proteomics_10percent_filtered_25_07_2025.rds"
)

predictor <- "30038_0_0"  # VO2max
formula <- "+ age_group + sex + ethnicity"

total_cols <- ncol(data_merged_filtered_proteomics)
end_col <- min(2973, total_cols)
if (total_cols < 51) stop("No hay suficientes columnas para seleccionar proteínas.")
proteins <- names(data_merged_filtered_proteomics)[51:end_col]

# ==========================
# 4. FORMAT VARIABLES
# ==========================
categorical_vars <- c("age_group", "sex", "ethnicity")
for (var in categorical_vars) {
  if (var %in% names(data_merged_filtered_proteomics)) {
    data_merged_filtered_proteomics[[var]] <- as.factor(data_merged_filtered_proteomics[[var]])
  } else {
    warning(paste("Variable categórica", var, "no encontrada en el dataset."))
  }
}

data_merged_filtered_proteomics[, proteins] <- lapply(
  data_merged_filtered_proteomics[, proteins, drop = FALSE],
  function(x) as.numeric(as.character(x))
)

# ==========================
# 5. RUN ANALYSIS
# ==========================
results_proteins <- recursive_lm(
  outcome = proteins,
  predictor = predictor,
  data = data_merged_filtered_proteomics,
  formula = formula,
  min_n = 30
)

# FDR
if (nrow(results_proteins) > 0) {
  results_proteins <- results_proteins %>%
    rename(p_value_raw = p.value) %>%
    mutate(
      p_value_fdr = p.adjust(p_value_raw, method = "fdr"),
      significant_fdr = p_value_fdr < 0.05
    )
}

# ==========================
# 6. SAVE MAIN RESULTS (.xlsx)
# ==========================
if (nrow(results_proteins) == 0) {
  cat("No models were generated. Check filters and function.\n")
} else {
  writexl::write_xlsx(
    results_proteins,
    "./upload/proteomics_associations.xlsx"
  )
}

# ==========================
# 7. SAMPLE SIZE SUMMARY (.xlsx)
# ==========================
covariates <- trimws(unlist(strsplit(formula, "\\+")))
covariates <- covariates[covariates != ""]
vars_for_max_N <- c(predictor, covariates)

max_N <- data_merged_filtered_proteomics %>%
  select(all_of(vars_for_max_N)) %>%
  filter(complete.cases(.)) %>%
  nrow()

if (nrow(results_proteins) > 0) {
  n_summary_proteins <- results_proteins %>%
    select(outcome, N) %>%
    distinct() %>%
    mutate(
      max_N = max_N,
      loss_abs = max_N - N,
      loss_rel_pct = round((max_N - N) / max_N * 100, 1)
    ) %>%
    arrange(desc(N))
  
  print(n_summary_proteins)
  writexl::write_xlsx(
    n_summary_proteins,
    "./upload/univariate_sample_size_per_protein.xlsx"
  )
}

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
