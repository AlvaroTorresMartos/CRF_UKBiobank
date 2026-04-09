

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  survival, broom, Hmisc, 
               ggplot2, tibble, ggtext, patchwork, forcats)

# Import data -----

covariates = read_rds("/mnt/project/data/processed/2025_07_25_covariates.RDS")

data = read_rds("/mnt/project/data/processed/validation_population_proteomic_signature.RDS")

data = data %>% 
  dplyr::filter(eid %in% covariates$eid) %>% 
  dplyr::select(eid, sex, z_scaled)

covariates = covariates %>% 
  dplyr::select(eid, 
                age_at_recruitment, starts_with("age_out"),
                age_5y_group, ethnicity, townsend_q5, education_group, 
                country, smoking_status, alc_group, fruit_veg_rec, satfat_3, bmi_group, 
                diab_fam, cvd_fam, htn_fam, cancer_fam, menopause, 
                medication_cholesterol, medication_bp, medication_diabetes, 
                medication_others, 
                t2dm, cvd, any_death, bladder, breast, colon, colorectal, 
                endometrial, lung, prostate)

data_def = data %>% 
  dplyr::left_join(covariates, by = "eid")


## T2DM ----
data_def = data_def %>% 
  dplyr::mutate(quantile_outcome = ntile(z_scaled, n = 5), 
                quantile_outcome = factor(quantile_outcome)
  )

surv_obj = with(data_def, Surv(time = age_at_recruitment,
                                   time2 = age_out_diab,
                                   event = t2dm == 1))

cox = coxph(surv_obj ~ quantile_outcome + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3   +
              diab_fam +
              medication_cholesterol + medication_bp + 
              strata(sex, bmi_group, menopause),
                          data = data_def)

cox.zph(cox)

summary(cox)

cox_t2dm = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_t2dm = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_t2dm = sum(t2dm == 1, na.rm = TRUE), 
    perc_t2dm = 100 * mean(t2dm == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_t2dm = cox_t2dm %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


refs = data.frame("term" = "quantile_outcome1", 
                  "estimate" = 1, 
                  "std.error" = 0, 
                  "statistic" = 0,
                  "p.value" = "", 
                  "conf.low" = 0, 
                  "conf.high" = 0)

cox_t2dm = rbind(refs, cox_t2dm)

## CVD ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_cvd,
                               event = cvd == 1))

# cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
#               townsend_q5 + education_group + country +
#               smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
#               cvd_fam + menopause +
#               medication_cholesterol + medication_bp + medication_others,
#                           data = data_def)

cox = coxph(surv_obj ~ quantile_outcome +
              townsend_q5 + country + alc_group + fruit_veg_rec + bmi_group +
               menopause +
              medication_cholesterol  + medication_bp + medication_others + 
              strata(sex, ethnicity, education_group, smoking_status, satfat_3, cvd_fam, medication_bp),
            data = data_def)



cox.zph(cox)

summary(cox)

cox_cvd = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_cvd = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_cvd = sum(cvd == 1, na.rm = TRUE), 
    perc_cvd = 100 * mean(cvd == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_cvd = cox_cvd %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_cvd = rbind(refs, cox_cvd)

## Death ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_death,
                               event = any_death == 1))

cox = coxph(surv_obj ~ quantile_outcome +  ethnicity +
               education_group + country +
               alc_group + fruit_veg_rec + satfat_3  +
              diab_fam + cvd_fam + cancer_fam +
              medication_cholesterol + medication_bp + medication_others + 
              strata(sex, townsend_q5, smoking_status, bmi_group, menopause),
            data = data_def)


cox.zph(cox)

summary(cox)

cox_death = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_death = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_death = sum(any_death == 1, na.rm = TRUE), 
    perc_death = 100 * mean(any_death == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_death = cox_death %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_death = rbind(refs, cox_death)

## Bladder ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_bladder,
                               event = bladder == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)


cox.zph(cox)

summary(cox)

cox_bladder = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_bladder = data_def %>%
  dplyr::group_by(quantile_outcome) %>%
  dplyr::summarise(
    count_total = n(),
    count_bladder = sum(bladder == 1, na.rm = TRUE),
    perc_bladder = 100 * mean(bladder == 1, na.rm = TRUE)
  ) %>%
  dplyr::select(-c(quantile_outcome))


cox_bladder = cox_bladder %>%
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3",
                            "quantile_outcome4", "quantile_outcome5"))


cox_bladder = rbind(refs, cox_bladder)


## Breast ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_breast,
                               event = breast == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)


cox.zph(cox)

summary(cox)

cox_breast = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_breast = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_breast = sum(breast == 1, na.rm = TRUE), 
    perc_breast = 100 * mean(breast == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_breast = cox_breast %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_breast = rbind(refs, cox_breast)


## Colon ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_colon,
                               event = colon == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)


cox.zph(cox)

summary(cox)

cox_colon = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_colon = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_colon = sum(colon == 1, na.rm = TRUE), 
    perc_colon = 100 * mean(colon == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_colon = cox_colon %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_colon = rbind(refs, cox_colon)

## Colorectal ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_colorectal,
                               event = colorectal == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)


cox.zph(cox)

summary(cox)

cox_colorectal = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_colorectal = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_colorectal = sum(colorectal == 1, na.rm = TRUE), 
    perc_colorectal = 100 * mean(colorectal == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_colorectal = cox_colorectal %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_colorectal = rbind(refs, cox_colorectal)

## Endometrial ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_endometrial,
                               event = endometrial == 1))

cox = coxph(surv_obj ~ quantile_outcome + ethnicity +
              townsend_q5 + education_group +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam +  medication_bp + medication_others + 
              strata(sex, menopause, medication_cholesterol),
            data = data_def)

cox.zph(cox)

summary(cox)

cox_endometrial = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_endometrial = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_endometrial = sum(endometrial == 1, na.rm = TRUE), 
    perc_endometrial = 100 * mean(endometrial == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_endometrial = cox_endometrial %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_endometrial = rbind(refs, cox_endometrial)

## Prostate ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_prostate,
                               event = prostate == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)



cox.zph(cox)

summary(cox)

cox_prostate = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_prostate = data_def %>%
  dplyr::group_by(quantile_outcome) %>%
  dplyr::summarise(
    count_total = n(),
    count_prostate = sum(prostate == 1, na.rm = TRUE),
    perc_prostate = 100 * mean(prostate == 1, na.rm = TRUE)
  ) %>%
  dplyr::select(-c(quantile_outcome))


cox_prostate = cox_prostate %>%
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3",
                            "quantile_outcome4", "quantile_outcome5"))


cox_prostate = rbind(refs, cox_prostate)

## Lung ----
surv_obj = with(data_def, Surv(time = age_at_recruitment,
                               time2 = age_out_lung,
                               event = lung == 1))

cox = coxph(surv_obj ~ quantile_outcome + sex + ethnicity +
              townsend_q5 + education_group + country +
              smoking_status + alc_group + fruit_veg_rec + satfat_3  + bmi_group +
              cancer_fam + menopause +
              medication_cholesterol + medication_bp + medication_others,
            data = data_def)

cox.zph(cox)

summary(cox)

cox_lung = broom::tidy(cox,  exponentiate = TRUE, conf.int = TRUE)

info_lung = data_def %>% 
  dplyr::group_by(quantile_outcome) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_lung = sum(lung == 1, na.rm = TRUE), 
    perc_lung = 100 * mean(lung == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_outcome))


cox_lung = cox_lung %>% 
  dplyr::filter(term %in% c("quantile_outcome2", "quantile_outcome3", 
                            "quantile_outcome4", "quantile_outcome5")) 


cox_lung = rbind(refs, cox_lung)

# Forest plots -----

# T2DM ----

results = cox_t2dm %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_t2dm)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_t2dm, perc_t2dm)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# CVD ----

results = cox_cvd %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_cvd)


forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_cvd, perc_cvd)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.3f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Bladder ----

results = cox_bladder %>%
  dplyr::mutate(HR = if_else(p.value == "", "Reference",
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score",
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3",
                                 term %in% "quantile_outcome4" ~ "Quintile 4",
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                ))
results = cbind(results, info_bladder)

forest_prote = results %>%
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"),
    axis.text.x = element_text(color = "black"))  +
  labs(x = "", y = "")

forest_prote





table1 = results %>%
  dplyr::mutate(term = paste0("**", term, "**"),
                term = forcats::fct_rev(term)) %>%
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(),
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**")

table1



table2 = results %>%
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(),
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**")

table2

table3 = results %>%
  dplyr::mutate(term = forcats::fct_rev(term),
                count_percentage = sprintf("%.0f (%.2f%%)", count_bladder, perc_bladder)) %>%
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(),
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**")

table3

table4 = results %>%
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(),
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**")

table4

table5 = results %>%
  dplyr::mutate(term = forcats::fct_rev(term),
                p.value = as.numeric(p.value),
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001",
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>%
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(),
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**")

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Any death ----

results = cox_death %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_death)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_death, perc_death)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Breast ----

results = cox_breast %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_breast)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_breast, perc_breast)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Colon ----

results = cox_colon %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_colon)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_colon, perc_colon)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Colorectal ----

results = cox_colorectal %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_colorectal)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_colorectal, perc_colorectal)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Endometrial ----

results = cox_endometrial %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_endometrial)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_endometrial, perc_endometrial)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Prostate ----

results = cox_prostate %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_prostate)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_prostate, perc_prostate)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

# Lung ----

results = cox_lung %>% 
  dplyr::mutate(HR = if_else(p.value == "", "Reference", 
                             sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
                term = case_when(term %in% "z_scaled" ~ "Continuous score", 
                                 term %in% "quantile_outcome1" ~ "Quintile 1",
                                 term %in% "quantile_outcome2" ~ "Quintile 2",
                                 term %in% "quantile_outcome3" ~ "Quintile 3", 
                                 term %in% "quantile_outcome4" ~ "Quintile 4", 
                                 term %in% "quantile_outcome5" ~ "Quintile 5"
                )) 
results = cbind(results, info_lung)

forest_prote = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 2, 4, 6)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomic signature**") 

table1



table2 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_lung, perc_lung)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)

