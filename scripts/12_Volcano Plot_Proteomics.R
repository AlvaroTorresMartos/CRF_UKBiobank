
library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)

data <- read_excel("/mnt/project/data/processed/proteomics_associations.xlsx")

data <- data %>%
  mutate(
    beta_coefficient = as.numeric(beta_coefficient),
    p_value_fdr = as.numeric(p_value_fdr)
  ) %>%
  filter(!is.na(beta_coefficient), !is.na(p_value_fdr), p_value_fdr > 0)

data <- data %>%
  mutate(
    neg_log10_p = -log10(p_value_fdr),
    significance_category = case_when(
      p_value_fdr < 0.05 & beta_coefficient > 0 ~ "Positive",
      p_value_fdr < 0.05 & beta_coefficient < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  )

data <- data %>%
  arrange(p_value_fdr) %>%
  mutate(custom_label = ifelse(p_value_fdr < 0.05 & row_number() <= 30, outcome, NA))

pdf("./upload/proteomics_volcano_plot.pdf", width = 8, height = 6)

ggplot(data, aes(x = beta_coefficient, y = neg_log10_p)) +
  geom_point(aes(color = significance_category), size = 2, alpha = 0.85) +
  scale_color_manual(values = c(
    "Positive" = "#D73027",
    "Negative" = "#4575B4",
    "Not significant" = "grey70"
  )) +
  scale_x_continuous(
    breaks = seq(-0.15, 0.15, by = 0.05),
    limits = c(-0.15, 0.15),
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression(-log[10]("FDR-adjusted P-value")),
    breaks = seq(0, 270, by = 50),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  geom_text_repel(
    data = subset(data, !is.na(custom_label)),
    aes(label = custom_label),
    size = 3,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.2,
    max.overlaps = 30
  ) +
  labs(
    x = "Beta coefficient",
    title = "Associations between VO₂max and proteins"
  ) +
  coord_cartesian(ylim = c(0, 270)) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
