#########################################################
# UNIVARIATE LINEAR REGRESSIONS FOR MULTIPLE METABOLITES
# Includes standardized and unstandardized betas
# Adjusted for confounders
# 17/02/2026
#########################################################

# ==========================
# 1. LOAD REQUIRED PACKAGES
# ==========================
install.packages("broom")
install.packages("dplyr")
install.packages("lm.beta")
install.packages("writexl")
install.packages("readxl")

library(broom)
library(dplyr)
library(lm.beta)
library(writexl)
library(readxl)

# ==========================
# 2. MAIN FUNCTION
# ==========================
recursive_lm <- function(outcome, predictor, data, formula, min_n = 50) {
  
  results <- lapply(outcome, function(y) {
    
    # Clean covariates from formula
    covariates <- trimws(unlist(strsplit(formula, "\\+")))
    covariates <- covariates[covariates != ""]
    
    # Define variables to keep in the model
    model_vars <- c(y, predictor, covariates)
    data_model <- data[complete.cases(data[, model_vars]), ]
    
    # Print N to debug
    cat("🔬 Metabolito:", y, "- N =", nrow(data_model), "\n")
    
    # Skip if too few observations
    if (nrow(data_model) < min_n) return(NULL)
    
    # Raw model
    object2 <- lm(
      as.formula(paste0(y, " ~ `", predictor, "`", formula)), data = data_model
    )
    
    # Extract raw coefficients + CI + N
    tidy_object2 <- broom::tidy(object2, conf.int = TRUE) %>%
      filter(term == predictor | term == paste0("`", predictor, "`"))
    
    # Check if predictor term was found
    if (nrow(tidy_object2) == 0) {
      cat("⚠️ No se encontró el término del predictor en:", y, "\n")
      return(NULL)
    }
    
    # Standardized model
    object1 <- lm.beta::lm.beta(
      lm(as.formula(paste0(y, " ~ `", predictor, "`", formula)), data = data_model)
    )
    
    # Extract standardized betas
    tidy_object1 <- broom::tidy(object1) %>%
      select(term, std_estimate) %>%
      rename(predictor = term, std_coefficient = std_estimate)
    
    # Complete tidy_object2 with additional columns
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
    
    # Combine results
    tidy_object <- merge(tidy_object1, tidy_object2, by = "predictor") %>%
      relocate(outcome, .before = predictor) %>%
      relocate(beta_coefficient, .before = std_coefficient)
    
    return(tidy_object)
  })
  
  # Bind all models
  final_result <- bind_rows(results)
  return(final_result)
}

# ==========================
# 3. LOAD DATA & DEFINE VARIABLES
# ==========================
data_merged_filtered_metabolomics <- readRDS(
  file.path("/discovery_cohort_metabolomics.rds")
)

predictor <- "30038_0_0"  # VO2max
formula <- "+ age_group + sex + ethnicity"

# Identify only metabolite columns starting with metabo_ and ending in _0_0
metabolites <- grep("^metabo_.*_0_0$", names(data_merged_filtered_metabolomics), value = TRUE)

# ==========================
# 3A. FORMAT VARIABLES
# ==========================
categorical_vars <- c("age_group", "sex", "ethnicity")
for (var in categorical_vars) {
  data_merged_filtered_metabolomics[[var]] <- as.factor(data_merged_filtered_metabolomics[[var]])
}

data_merged_filtered_metabolomics[, metabolites] <- lapply(
  data_merged_filtered_metabolomics[, metabolites],
  function(x) as.numeric(as.character(x))
)

# ==========================
# 4. RUN ANALYSIS
# ==========================
results <- recursive_lm(
  outcome = metabolites,
  predictor = predictor,
  data = data_merged_filtered_metabolomics,
  formula = formula,
  min_n = 50
)

# FDR adjust
if (!is.null(results) && nrow(results) > 0) {
  results <- results %>%
    rename(p_value_raw = p.value) %>%
    mutate(
      p_value_fdr = p.adjust(p_value_raw, method = "fdr"),
      significant_fdr = p_value_fdr < 0.05
    )
}

# ==========================
# 4B. ADD METABOLITE NAMES (title) BEFORE EXPORT
# ==========================
metabo_dict <- readxl::read_xlsx(
  "/metabo_names.xlsx"
) %>%
  dplyr::select(field_id, title, Group, Subgroup, variable_name) %>%
  dplyr::mutate(variable_name = trimws(variable_name))

if (!is.null(results) && nrow(results) > 0) {
  results <- results %>%
    dplyr::mutate(outcome = trimws(outcome)) %>%
    dplyr::left_join(metabo_dict, by = c("outcome" = "variable_name")) %>%
    dplyr::rename(outcome_title = title) %>%
    dplyr::relocate(outcome_title, Group, Subgroup, .after = outcome)
}

# ==========================
# 5. SAVE MAIN RESULTS (.xlsx)
# ==========================
if (is.null(results) || nrow(results) == 0) {
  cat("⚠️ No se generó ningún modelo. Revisa los filtros y la función.\n")
} else {
  writexl::write_xlsx(results, "/metabolomics_associations.xlsx")
}

# ==========================
# 6. SAMPLE SIZE SUMMARY (.xlsx)
# ==========================
covariates <- trimws(unlist(strsplit(formula, "\\+")))
covariates <- covariates[covariates != ""]
vars_for_max_N <- c(predictor, covariates)

max_N <- data_merged_filtered_metabolomics %>%
  select(all_of(vars_for_max_N)) %>%
  filter(complete.cases(.)) %>%
  nrow()

if (!is.null(results) && nrow(results) > 0) {
  n_summary <- results %>%
    select(outcome, N) %>%
    distinct() %>%
    mutate(
      max_N = max_N,
      loss_abs = max_N - N,
      loss_rel_pct = round((max_N - N) / max_N * 100, 1)
    ) %>%
    arrange(desc(N))
  
  print(n_summary)
  writexl::write_xlsx(n_summary, "/metabolomics_associations_summary.xlsx")
}
