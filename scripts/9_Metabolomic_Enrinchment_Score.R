# =========================================================
# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026
# =========================================================

library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(readxl)

# ============================
# 1) Import linear regression results
# ============================
assoc <- read_excel(
  "/metabolomics_associations.xlsx",
  sheet = "Sheet1"
)

# ============================
# 2) Create ranking for fgsea (USE outcome_title + beta_coefficient)
#     - Convert comma decimals
#     - Ensure unique names (mean if duplicated)
# ============================
assoc2 <- assoc %>%
  mutate(
    outcome_title = trimws(as.character(outcome_title)),
    beta_coefficient = as.numeric(gsub(",", ".", as.character(beta_coefficient)))
  ) %>%
  filter(!is.na(outcome_title), outcome_title != "", !is.na(beta_coefficient))

ranking <- assoc2 %>%
  group_by(outcome_title) %>%
  summarise(beta = mean(beta_coefficient), .groups = "drop") %>%
  select(outcome_title, beta) %>%
  tibble::deframe()

ranking <- sort(ranking, decreasing = TRUE)

stopifnot(is.numeric(ranking))
stopifnot(!any(is.na(ranking)))
stopifnot(!any(duplicated(names(ranking))))

# ============================
# 3) Define metabolite functional classes (DESCRIPTIVE NAMES)
#    Based on your real outcome_title list (no underscores)
# ============================

met_classes <- list(
  
  Cholesterol = c(
    "Total cholesterol",
    "Total cholesterol minus HDL-C",
    "Remnant cholesterol (non-HDL, non-LDL -cholesterol)",
    "VLDL cholesterol",
    "Cholesterol in IDL",
    "Clinical LDL cholesterol",
    "LDL cholesterol",
    "HDL cholesterol",
    
    "Cholesterol in chylomicrons and extremely large VLDL",
    "Cholesterol in very large VLDL",
    "Cholesterol in large VLDL",
    "Cholesterol in medium VLDL",
    "Cholesterol in small VLDL",
    "Cholesterol in very small VLDL",
    
    "Cholesterol in large LDL",
    "Cholesterol in medium LDL",
    "Cholesterol in small LDL",
    
    "Cholesterol in very large HDL",
    "Cholesterol in large HDL",
    "Cholesterol in medium HDL",
    "Cholesterol in small HDL",
    
    "Cholesterol to total lipids ratio in chylomicrons and extremely large VLDL",
    "Cholesterol to total lipids ratio in very large VLDL",
    "Cholesterol to total lipids ratio in large VLDL",
    "Cholesterol to total lipids ratio in medium VLDL",
    "Cholesterol to total lipids ratio in small VLDL",
    "Cholesterol to total lipids ratio in very small VLDL",
    "Cholesterol to total lipids ratio in IDL",
    "Cholesterol to total lipids ratio in large LDL",
    "Cholesterol to total lipids ratio in medium LDL",
    "Cholesterol to total lipids ratio in small LDL",
    "Cholesterol to total lipids ratio in very large HDL",
    "Cholesterol to total lipids ratio in large HDL",
    "Cholesterol to total lipids ratio in medium HDL",
    "Cholesterol to total lipids ratio in small HDL"
  ),
  
  Triglycerides = c(
    "Total triglycerides",
    "Triglycerides in VLDL",
    "Triglycerides in IDL",
    "Triglycerides in LDL",
    "Triglycerides in HDL",
    
    "Triglycerides in chylomicrons and extremely large VLDL",
    "Triglycerides in very large VLDL",
    "Triglycerides in large VLDL",
    "Triglycerides in medium VLDL",
    "Triglycerides in small VLDL",
    "Triglycerides in very small VLDL",
    
    "Triglycerides in large LDL",
    "Triglycerides in medium LDL",
    "Triglycerides in small LDL",
    
    "Triglycerides in very large HDL",
    "Triglycerides in large HDL",
    "Triglycerides in medium HDL",
    "Triglycerides in small HDL",
    
    "Triglycerides to total lipids ratio in chylomicrons and extremely large VLDL",
    "Triglycerides to total lipids ratio in very large VLDL",
    "Triglycerides to total lipids ratio in large VLDL",
    "Triglycerides to total lipids ratio in medium VLDL",
    "Triglycerides to total lipids ratio in small VLDL",
    "Triglycerides to total lipids ratio in very small VLDL",
    "Triglycerides to total lipids ratio in IDL",
    "Triglycerides to total lipids ratio in large LDL",
    "Triglycerides to total lipids ratio in medium LDL",
    "Triglycerides to total lipids ratio in small LDL",
    "Triglycerides to total lipids ratio in very large HDL",
    "Triglycerides to total lipids ratio in large HDL",
    "Triglycerides to total lipids ratio in medium HDL",
    "Triglycerides to total lipids ratio in small HDL",
    
    "Ratio of triglycerides to phosphoglycerides"
  ),
  
  Phospholipids = c(
    "Total phospholipids in lipoprotein particles",
    
    "Phospholipids in VLDL",
    "Phospholipids in IDL",
    "Phospholipids in LDL",
    "Phospholipids in HDL",
    
    "Phospholipids in chylomicrons and extremely large VLDL",
    "Phospholipids in very large VLDL",
    "Phospholipids in large VLDL",
    "Phospholipids in medium VLDL",
    "Phospholipids in small VLDL",
    "Phospholipids in very small VLDL",
    
    "Phospholipids in large LDL",
    "Phospholipids in medium LDL",
    "Phospholipids in small LDL",
    
    "Phospholipids in very large HDL",
    "Phospholipids in large HDL",
    "Phospholipids in medium HDL",
    "Phospholipids in small HDL",
    
    "Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL",
    "Phospholipids to total lipids ratio in very large VLDL",
    "Phospholipids to total lipids ratio in large VLDL",
    "Phospholipids to total lipids ratio in medium VLDL",
    "Phospholipids to total lipids ratio in small VLDL",
    "Phospholipids to total lipids ratio in very small VLDL",
    "Phospholipids to total lipids ratio in IDL",
    "Phospholipids to total lipids ratio in large LDL",
    "Phospholipids to total lipids ratio in medium LDL",
    "Phospholipids to total lipids ratio in small LDL",
    "Phospholipids to total lipids ratio in very large HDL",
    "Phospholipids to total lipids ratio in large HDL",
    "Phospholipids to total lipids ratio in medium HDL",
    "Phospholipids to total lipids ratio in small HDL",
    
    "Phosphatidylcholines",
    "Phosphoglycerides",
    "Sphingomyelins",
    "Total cholines"
  ),
  
  Cholesteryl_esters = c(
    "Total esterified cholesterol",
    
    "Cholesteryl esters in VLDL",
    "Cholesteryl esters in IDL",
    "Cholesteryl esters in LDL",
    "Cholesteryl esters in HDL",
    
    "Cholesteryl esters in chylomicrons and extremely large VLDL",
    "Cholesteryl esters in very large VLDL",
    "Cholesteryl esters in large VLDL",
    "Cholesteryl esters in medium VLDL",
    "Cholesteryl esters in small VLDL",
    "Cholesteryl esters in very small VLDL",
    
    "Cholesteryl esters in large LDL",
    "Cholesteryl esters in medium LDL",
    "Cholesteryl esters in small LDL",
    
    "Cholesteryl esters in very large HDL",
    "Cholesteryl esters in large HDL",
    "Cholesteryl esters in medium HDL",
    "Cholesteryl esters in small HDL",
    
    "Cholesteryl esters to total lipids ratio in chylomicrons and extremely large VLDL",
    "Cholesteryl esters to total lipids ratio in very large VLDL",
    "Cholesteryl esters to total lipids ratio in large VLDL",
    "Cholesteryl esters to total lipids ratio in medium VLDL",
    "Cholesteryl esters to total lipids ratio in small VLDL",
    "Cholesteryl esters to total lipids ratio in very small VLDL",
    "Cholesteryl esters to total lipids ratio in IDL",
    "Cholesteryl esters to total lipids ratio in large LDL",
    "Cholesteryl esters to total lipids ratio in medium LDL",
    "Cholesteryl esters to total lipids ratio in small LDL",
    "Cholesteryl esters to total lipids ratio in very large HDL",
    "Cholesteryl esters to total lipids ratio in large HDL",
    "Cholesteryl esters to total lipids ratio in medium HDL",
    "Cholesteryl esters to total lipids ratio in small HDL"
  ),
  
  Free_cholesterol = c(
    "Total free cholesterol",
    
    "Free cholesterol in VLDL",
    "Free cholesterol in IDL",
    "Free cholesterol in LDL",
    "Free cholesterol in HDL",
    
    "Free cholesterol in chylomicrons and extremely large VLDL",
    "Free cholesterol in very large VLDL",
    "Free cholesterol in large VLDL",
    "Free cholesterol in medium VLDL",
    "Free cholesterol in small VLDL",
    "Free cholesterol in very small VLDL",
    
    "Free cholesterol in large LDL",
    "Free cholesterol in medium LDL",
    "Free cholesterol in small LDL",
    
    "Free cholesterol in very large HDL",
    "Free cholesterol in large HDL",
    "Free cholesterol in medium HDL",
    "Free cholesterol in small HDL",
    
    "Free cholesterol to total lipids ratio in chylomicrons and extremely large VLDL",
    "Free cholesterol to total lipids ratio in very large VLDL",
    "Free cholesterol to total lipids ratio in large VLDL",
    "Free cholesterol to total lipids ratio in medium VLDL",
    "Free cholesterol to total lipids ratio in small VLDL",
    "Free cholesterol to total lipids ratio in very small VLDL",
    "Free cholesterol to total lipids ratio in IDL",
    "Free cholesterol to total lipids ratio in large LDL",
    "Free cholesterol to total lipids ratio in medium LDL",
    "Free cholesterol to total lipids ratio in small LDL",
    "Free cholesterol to total lipids ratio in very large HDL",
    "Free cholesterol to total lipids ratio in large HDL",
    "Free cholesterol to total lipids ratio in medium HDL",
    "Free cholesterol to total lipids ratio in small HDL"
  ),
  
  Total_lipids = c(
    "Total lipids in lipoprotein particles",
    
    "Total lipids in VLDL",
    "Total lipids in IDL",
    "Total lipids in LDL",
    "Total lipids in HDL",
    
    "Total lipids in chylomicrons and extremely large VLDL",
    "Total lipids in very large VLDL",
    "Total lipids in large VLDL",
    "Total lipids in medium VLDL",
    "Total lipids in small VLDL",
    "Total lipids in very small VLDL",
    
    "Total lipids in large LDL",
    "Total lipids in medium LDL",
    "Total lipids in small LDL",
    
    "Total lipids in very large HDL",
    "Total lipids in large HDL",
    "Total lipids in medium HDL",
    "Total lipids in small HDL"
  ),
  
  LP_concentration = c(
    "Total concentration of lipoprotein particles",
    "Concentration of chylomicrons and extremely large VLDL particles",
    
    "Concentration of VLDL particles",
    "Concentration of IDL particles",
    "Concentration of LDL particles",
    "Concentration of HDL particles",
    
    "Concentration of very large VLDL particles",
    "Concentration of large VLDL particles",
    "Concentration of medium VLDL particles",
    "Concentration of small VLDL particles",
    "Concentration of very small VLDL particles",
    
    "Concentration of large LDL particles",
    "Concentration of medium LDL particles",
    "Concentration of small LDL particles",
    
    "Concentration of very large HDL particles",
    "Concentration of large HDL particles",
    "Concentration of medium HDL particles",
    "Concentration of small HDL particles",
    
    "Concentration of very large VLDL particles",
    "Concentration of very large HDL particles"
  ),
  
  ApoLP_size = c(
    "Apolipoprotein A1",
    "Apolipoprotein B",
    "Ratio of apolipoprotein B to apolipoprotein A1",
    "Average diameter for VLDL particles",
    "Average diameter for LDL particles",
    "Average diameter for HDL particles"
  ),
  
  Fatty_acids = c(
    "Total fatty acids",
    "Degree of unsaturation",
    "Omega-3 fatty acids",
    "Omega-6 fatty acids",
    "Polyunsaturated fatty acids",
    "Monounsaturated fatty acids",
    "Saturated fatty acids",
    "Linoleic acid",
    "Docosahexaenoic acid",
    
    "Ratio of omega-3 fatty acids to total fatty acids",
    "Ratio of omega-6 fatty acids to total fatty acids",
    "Ratio of omega-6 fatty acids to omega-3 fatty acids",
    "Ratio of polyunsaturated fatty acids to total fatty acids",
    "Ratio of monounsaturated fatty acids to total fatty acids",
    "Ratio of saturated fatty acids to total fatty acids",
    "Ratio of linoleic acid to total fatty acids",
    "Ratio of docosahexaenoic acid to total fatty acids",
    "Ratio of polyunsaturated fatty acids to monounsaturated fatty acids"
  ),
  
  Amino_acids = c(
    "Alanine",
    "Spectrometer-corrected alanine",
    "Glutamine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Leucine",
    "Valine",
    "Phenylalanine",
    "Tyrosine",
    "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)"
  ),
  
  Glycolysis = c(
    "Glucose",
    "Lactate",
    "Pyruvate",
    "Citrate",
    "Glucose-lactate"
  ),
  
  Ketone_bodies = c(
    "3-Hydroxybutyrate",
    "Acetate",
    "Acetoacetate",
    "Acetone"
  ),
  
  Fluid_balance = c(
    "Creatinine",
    "Albumin"
  ),
  
  Inflammation = c(
    "Glycoprotein acetyls"
  )
)

# ============================
# 3b) Safety checks: all pathway terms must exist in ranking
# ============================
all_terms <- unique(unlist(met_classes))
overlap <- sum(all_terms %in% names(ranking))
cat("Overlap terms:", overlap, "of", length(all_terms),
    sprintf("(%.1f%%)\n", 100*overlap/length(all_terms)))

missing <- setdiff(all_terms, names(ranking))
if (length(missing) > 0) {
  cat("Missing terms (first 50):\n")
  print(head(missing, 50))
  stop("Fix met_classes terms or outcome_title naming before running fgsea.")
}

# ============================
# 4) Execute fgsea (multilevel)
#    Note: setting maxSize=50 may result in 0 pathways passing filters.
# ============================
fgsea_results <- fgseaMultilevel(
  pathways = met_classes,
  stats    = ranking,
  minSize  = 1,
  maxSize  = 300
)

# Sort by adjusted significance
fgsea_results <- fgsea_results[order(fgsea_results$padj), ]

# Convert leadingEdge into a single string for readability
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))

# Prettify pathway names
fgsea_results <- fgsea_results %>%
  mutate(pathway = gsub("_", " ", pathway)) %>%
  mutate(Significance = ifelse(padj < 0.05, "*", ""))

# ============================
# 5) Plot (guard-rail to avoid errors if there are 0 rows)
# ============================
if (nrow(fgsea_results) == 0) {
  msea_plot <- ggplot() +
    theme_void() +
    ggtitle("No fgsea results (0 pathways passed minSize/maxSize filters)")
} else {
  msea_plot <- fgsea_results %>%
    mutate(direction = if_else(NES > 0, "Up-Regulated", "Down-Regulated")) %>%
    ggplot(aes(
      y = reorder(pathway, NES),
      x = NES,
      color = padj,   
      size = size
    )) +
    facet_grid(. ~ direction, scales = "free_x", space = "free_x") +
    scale_color_gradientn(name = "Adj. P-value", colors = c("red", "purple", "blue")) +
    scale_size(name = "Count", range = c(5, 12)) +
    geom_point() +
    labs(
      title = "Metabolite Set Enrichment Analysis (LM Betas)",
      y = "",
      x = "Enrichment score (NES)"
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "cm"),
      strip.text = element_text(size = 20),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

msea_plot

# ============================
# 6) Save outputs (PDF + Excel)
# ============================

# Output folder (your path)
out_dir <- "/results_new"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Save as PDF
ggsave(
  filename = file.path(out_dir, "/Metabolomics_Enrichment_Betas.pdf"),
  plot = msea_plot,
  width = 10,
  height = 8,
  device = cairo_pdf  # <- comment this line out if Cairo is not available
)

# ---- Export results to Excel
library(writexl)

write_xlsx(
  fgsea_results,
  path = file.path(out_dir, "/Metabolomics_Enrichment_Betas.xlsx")
)
