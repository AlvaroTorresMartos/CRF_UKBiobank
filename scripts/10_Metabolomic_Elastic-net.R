# =========================================================
# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026
# =========================================================

# Load packages ------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  readr, dplyr, janitor, stringr, naniar, broom,
  tidyr, tibble, Hmisc, glmnet, descr, fastDummies, ggplot2,
  caret, MLmetrics, foreign, parallel, doParallel,
  imputeLCMD, magrittr, foreach
)

# Import data -----
t1 <- Sys.time()

metabolomics <- read_rds("/mnt/project/data/processed/Biomarkers_all.rds")

covariates <- read_rds("/mnt/project/data/processed/2025_07_25_covariates.RDS")

ids <- read_rds("/mnt/project/data/processed/discovery_cohort_metabolomics.rds")

names <- read_tsv(
  "/mnt/project/Showcase metadata/fieldsum.tsv",
  show_col_types = FALSE
)

t2 <- Sys.time()
t2 - t1
rm(t1, t2)

# Discovery population ----
ids <- ids %>% 
  filter(eid %in% covariates$eid) %>% 
  select(eid, sex, `30038_0_0`) %>% 
  rename(vo2max = `30038_0_0`)

id_correct <- ids$eid

metabolomics <- metabolomics %>% 
  filter(eid %in% id_correct)

# Preprocessing ----
# Adapt dictionary to match format "20280-0.0"

names <- names %>%
  mutate(
    field_id = paste0(field_id, "-0.0"),
    title    = make_clean_names(title)
  )

metabolites <- names$field_id[c(1176:1177, 2313:2561)]
metabolites_label <- names$title[c(1176:1177, 2313:2561)]

# Select and rename
metabolomics <- metabolomics %>%
  select(eid, all_of(metabolites)) %>%
  rename_with(~ c("eid", metabolites_label))

# Missingness
missing_values_case <- naniar::miss_case_summary(metabolomics)
missing_values_vars <- naniar::miss_var_summary(metabolomics)

# Imputation
qrilc_object <- imputeLCMD::impute.QRILC(dataSet.mvs = metabolomics)
metabolomics_imp <- qrilc_object[[1]]

# Join datasets
covariates <- covariates %>% 
  select(eid, age_at_recruitment, ethnicity)

data <- ids %>% 
  left_join(covariates, by = "eid")

# Missingness in covariates
missing_values_case <- naniar::miss_case_summary(data)
missing_values_vars <- naniar::miss_var_summary(data)

# Dummy variables
data <- fastDummies::dummy_columns(
  data,
  select_columns     = "ethnicity",
  remove_first_dummy = TRUE
) %>% 
  select(-ethnicity)

# Merge metabolomics
data <- data %>% 
  left_join(metabolomics_imp, by = "eid") %>% 
  relocate(vo2max, .after = last_col())

# Scale metabolites
data <- data %>% 
  mutate(across(all_of(metabolites_label), scale))

sum(is.na(data))

# [1] 0

##########################################################
#saveRDS(data, "/metabolomics_dataset_elasticnet_crf.RDS")
##########################################################

# Training Predictive modelling (elasticnet) -----

data = readr::read_rds("/mnt/project/data/processed/metabolomics_dataset_elasticnet_crf.RDS")

data$age_at_recruitment = scale(data$age_at_recruitment)
data$vo2max = as.numeric(scale(data$vo2max))

## Set-ups -----
set.seed(123456789)
folds = caret::createFolds(data$vo2max, k = 10,
                        list = TRUE, returnTrain =  TRUE)

## Experimental design ------
set.seed(123456789)
control = trainControl(method = "cv", 
                       number = 10,
                       savePredictions = "all",
                       index = folds,
                       allowParallel = TRUE
)

### Parallel process and run the elasticnet model 
parallel::detectCores()
# [1] 22
library(doParallel)
getDoParWorkers()
library(doParallel)
cl <- makePSOCKcluster(65)
registerDoParallel(cl)
# https://stackoverflow.com/questions/41117127/r-parallel-programming-error-in-task-1-failed-could-not-find-function
clusterCall(cl, function() library(magrittr))
t1 = Sys.time()
set.seed(123456789)
elasticnet = train(x = data[, c(2:258)], y = data$vo2max, 
                         method = "glmnet",   
                         penalty.factor = c(rep(0, 6), rep(1, 251)),
                         metric = "RMSE",
                         maximize =  FALSE,
                         trControl = control, 
                         tuneGrid = expand.grid(alpha = 0.5, #seq(0, 1, by = 0.1), 
                                                lambda = seq(0, 1, by = 0.1))
                         
)
registerDoSEQ()
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1
# Time difference of 12.3476 secs
results = elasticnet$results


###########################################################
## Save the model as R object ----
#saveRDS(elasticnet, "/elasticnet_metabolomics2.RDS")
############################################################

# Final Predictive Model ------
elasticnet = readr::read_rds("/mnt/project/data/processed/elasticnet_metabolomics2.RDS")

set.seed(123456789)
control = trainControl(method = "none")

t1 = Sys.time()
set.seed(123456789)
elasticnet_final = train(x = data[, c(2:258)], y = data$vo2max, 
                         method = "glmnet",
                         penalty.factor = c(rep(0, 6), rep(1, 251)),
                         trControl = control,
                         tuneGrid = elasticnet$bestTune
)
t2 = Sys.time()
t2 - t1

## Save the model as R object ----
#saveRDS(elasticnet_final, "/final_elasticnet_metabolomics2.RDS")

elasticnet_final <- read_rds("/mnt/project/data/processed/final_elasticnet_metabolomics2.RDS")

pred = predict(elasticnet_final, data[, 2:258])

RMSE(pred, data$vo2max)
# 0.7062737
MAE(pred, data$vo2max)
# 0.5319059
R2 <- cor(pred, data$vo2max)^2
R2
# 0.5011707

coef_elasticnet = as.matrix(coef(elasticnet_final$finalModel, elasticnet$bestTune$lambda))
coef_elasticnet = coef_elasticnet %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "Coefficient")  


covariates = data %>% dplyr::select(sex:ethnicity_9) %>%
  colnames()

coef_elasticnet = coef_elasticnet %>%
  slice(which(coef_elasticnet$Coefficient %in% c(covariates, "(Intercept)")  == FALSE))  %>%  
  filter(s0 != 0) 


#saveRDS(coef_elasticnet, "/coeficientes_huella_metabolomics.RDS")
coef_elasticnet <-readRDS("/mnt/project/data/processed/coeficientes_huella_metabolomics.RDS")
# Coefficient plot -----
coef_elasticnet %>% 
  dplyr::rename(value = s0) %>% 
  dplyr::mutate(abs_value = abs(value)) %>%
  dplyr::slice_max(order_by = abs_value, n = 30) %>%
  dplyr::mutate(direction = ifelse(value > 0, "+", "-")) %>% 
  ggplot(aes(x = value, 
             y = reorder(Coefficient, value), 
             fill = direction)) + 
  geom_col() + 
  labs(x = "Regularized regression coefficients", y = "") + 
  scale_fill_manual(values = c("#A7C7E7", "#F7A1A1")) +
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 14))



colnames(coef_elasticnet)

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
