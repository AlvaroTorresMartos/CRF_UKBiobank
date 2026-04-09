# =========================================================
# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026
# =========================================================

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, janitor, stringr, naniar, broom, 
               tidyr, tibble, Hmisc, glmnet, descr, fastDummies, ggplot2, 
               caret, MLmetrics, foreign, parallel, doParallel, readr, purrr)

# Import data -----
t1 = Sys.time()
metabolomics = read_rds("/mnt/project/data/processed/Biomarkers_all.rds")

covariates = read_rds("/mnt/project/data/processed/2025_07_25_covariates.RDS")

ids_disc = read_rds("/mnt/project/data/processed//discovery_cohort_metabolomics.rds")

ids_vali = read_rds("/mnt/project/data/processed//validation_cohort_metabolomics.rds")

names = readr::read_tsv("/mnt/project/Showcase metadata/fieldsum.tsv", show_col_types = FALSE)


t2 = Sys.time()
t2-t1
# Time difference of 1.794389 mins
rm(t1, t2)

# Discovery and validation populations ----
nrow(ids_disc)
# 54321
ids_disc = ids_disc %>% 
  dplyr::filter(eid %in% covariates$eid) %>% 
  dplyr::select(eid, sex, `30038_0_0`, ) %>% 
  dplyr::rename(vo2max = `30038_0_0`)
nrow(ids_disc)
# 54321

nrow(ids_vali)
# 354222
ids_vali = ids_vali %>% 
  dplyr::filter(eid %in% covariates$eid) %>% 
  dplyr::select(eid, sex, `30038_0_0`, ) %>% 
  dplyr::rename(vo2max = `30038_0_0`)
nrow(ids_vali)
# 354222


id_correct1 = ids_disc$eid
id_correct2 = ids_vali$eid

nrow(metabolomics)
# 502134
metabolomics1 = metabolomics %>% 
  dplyr::filter(eid %in% id_correct1)
nrow(metabolomics1)
# 54321

nrow(metabolomics)
# 502134
metabolomics2 = metabolomics %>% 
  dplyr::filter(eid %in% id_correct2)
nrow(metabolomics2)
# 354222

# Preprocessing -----
## Change the format of column names ----
metabolomics1 <- metabolomics1 %>%
  dplyr::rename_with(~ .x %>%
                       stringr::str_replace("^X", "") %>%        
                       stringr::str_replace_all("\\.", "-") %>%  
                       stringr::str_replace("-0-0$", "-0.0"))   

metabolomics2 <- metabolomics2 %>%
  dplyr::rename_with(~ .x %>%
                       stringr::str_replace("^X", "") %>%
                       stringr::str_replace_all("\\.", "-") %>%
                       stringr::str_replace("-0-0$", "-0.0"))

#Adapt dictionary to match format "20280-0.0"
names <- names %>%
  dplyr::mutate(
    field_id = paste0(field_id, "-0.0"),
    title    = janitor::make_clean_names(title)
  )

metabolites <- names$field_id[c(1176:1177, 2313:2561)]
metabolites_label <- names$title[c(1176:1177, 2313:2561)]
metabolites_label = c("eid", metabolites_label)

## Select the id and metabolites and rename -----
metabolomics1 = metabolomics1 %>%
  dplyr::select(eid, all_of(metabolites)) %>%
  dplyr::rename_with(~ metabolites_label) 

metabolomics2 = metabolomics2 %>%
  dplyr::select(eid, all_of(metabolites)) %>%
  dplyr::rename_with(~ metabolites_label)

## Evaluate the missingness in metabolomics data -----
missing_values_case = naniar::miss_case_summary(metabolomics1)
missing_values_vars = naniar::miss_var_summary(metabolomics1)

missing_values_case2 = naniar::miss_case_summary(metabolomics2)
missing_values_vars2 = naniar::miss_var_summary(metabolomics2)

case_with_acceptable_nas = missing_values_case2 %>% 
  dplyr::filter(pct_miss < 10)

nrow(metabolomics2)
#  354222
metabolomics2 = metabolomics2 %>% 
  dplyr::slice(case_with_acceptable_nas$case)
nrow(metabolomics2)
#  353908

missing_values_case2 = naniar::miss_case_summary(metabolomics2)
missing_values_vars2 = naniar::miss_var_summary(metabolomics2)

## Imputation -----
qrilc_object = imputeLCMD::impute.QRILC(dataSet.mvs = metabolomics1)
metabolomics_imp1 = qrilc_object[[1]] 

qrilc_object = imputeLCMD::impute.QRILC(dataSet.mvs = metabolomics2)
metabolomics_imp2 = qrilc_object[[1]] 

## Joining the datasets -----
covariates = covariates %>% 
  dplyr::select(eid, age_at_recruitment, ethnicity)

data_disc = ids_disc %>% 
  dplyr::left_join(covariates, by = "eid")

data_vali = ids_vali %>% 
  dplyr::left_join(covariates, by = "eid")

## Evaluate the missingness in covariates data  -----
missing_values_case = naniar::miss_case_summary(data_disc)
missing_values_vars = naniar::miss_var_summary(data_disc)

missing_values_case = naniar::miss_case_summary(data_vali)
missing_values_vars = naniar::miss_var_summary(data_vali)

## Dummies variables -----
data_disc = fastDummies::dummy_columns(
  data_disc,
  select_columns     = "ethnicity",
  remove_first_dummy = TRUE) %>% 
  dplyr::select(-c(ethnicity))

data_vali = fastDummies::dummy_columns(
  data_vali,
  select_columns     = "ethnicity",
  remove_first_dummy = TRUE) %>% 
  dplyr::select(-c(ethnicity))

## Scaling the metabolomic variables -----
metabolites_label = metabolites_label[2:length(metabolites_label)]

metabolomics_imp1 = metabolomics_imp1 %>% 
  dplyr::mutate(across(all_of(metabolites_label), scale))

metabolomics_imp2 = metabolomics_imp2 %>% 
  dplyr::mutate(across(all_of(metabolites_label), scale))

data_disc = data_disc %>% 
  dplyr::left_join(metabolomics_imp1, by = "eid") %>% 
  dplyr::relocate(vo2max, .after = last_col())

data_vali = data_vali %>% 
  dplyr::right_join(metabolomics_imp2, by = "eid") %>% 
  dplyr::select(-c(vo2max))

data_disc$age_at_recruitment = scale(data_disc$age_at_recruitment)
data_vali$age_at_recruitment = scale(data_vali$age_at_recruitment)

data_disc$vo2max = as.numeric(scale(data_disc$vo2max))

sum(is.na(data_disc))
# [1] 0

sum(is.na(data_vali))
# [1] 0

# Final model elastinet -----
elasticnet_final = readr::read_rds("/mnt/project/data/processed/final_elasticnet_metabolomics2.RDS")

metabolomic_signature_disc = predict(elasticnet_final, data_disc[, c(2:258)])
metabolomic_signature_vali = predict(elasticnet_final, data_vali[, c(2:258)])

mean = mean(metabolomic_signature_disc)
sd = sd(metabolomic_signature_disc)

data_disc$met_signature = metabolomic_signature_disc
data_vali$met_signature = metabolomic_signature_vali

data_disc$met_signature %>% hist()
data_vali$met_signature %>% hist()

data_disc = data_disc %>% 
  dplyr::mutate(z_scaled = met_signature)

data_vali = data_vali %>% 
  dplyr::mutate(z_scaled = met_signature)

# data_disc = data_disc %>% 
#   dplyr::mutate(z_scaled = scale(met_signature))
# 
# data_vali = data_vali %>% 
#   dplyr::mutate(z_scaled = (met_signature-mean)/sd)

data_disc$z_scaled %>% hist()

data_vali$z_scaled %>% hist()

# disc_popu = disc_popu %>% 
#   dplyr::select(eid, z_scaled)
# 
# vali_popu = vali_popu %>% 
#   dplyr::select(eid, z_scaled)

saveRDS(data_disc, "./upload/discovery_metabolomics_signature.RDS")
saveRDS(data_vali, "./upload/validation_metabolomics_signature.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")