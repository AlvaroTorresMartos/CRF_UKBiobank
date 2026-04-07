# =========================================================
# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026
# =========================================================

library(dplyr)
library(readr)
library(ggplot2)

path_data_metab  <- "/metabolomics_dataset_elasticnet_crf.RDS"
path_model_metab <- "/final_elasticnet_metabolomics2.RDS"

data  <- read_rds(path_data_metab)
model <- read_rds(path_model_metab)

# CRF scaled
data$vo2max <- as.numeric(scale(as.numeric(data$vo2max)))

# Signature (z_scaled)
pred <- as.numeric(predict(model, newdata = data[, 2:258]))

metabolomics <- tibble(
  vo2max   = data$vo2max,
  z_scaled = as.numeric(scale(pred))
)

# =========================================================
# PLOT
# =========================================================

metabolomics %>% 
  dplyr::rename(CRF = vo2max) %>% 
  ggplot(aes(x = z_scaled, y = CRF)) + 
  geom_point()  + 
  geom_smooth(method = "lm") + 
  labs(x = "Metabolomic CRF signature", y = "CRF") + 
  ggtitle(paste("R²:", 0.50 ,
                " | RMSE:", 0.71 ,
                " | MAE:", 0.53)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18))








