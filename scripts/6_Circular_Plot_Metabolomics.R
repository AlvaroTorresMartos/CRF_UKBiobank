# Project: UK Biobank and cardiorespiratory fitness
# Date: 17/02/2026

find("coord_radial")

# =========================================================
# 0) Load packages ----
# =========================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, readxl, janitor, stringr, tibble, tidyr, purrr, scales,
  ggplot2, Cairo
)

# =========================================================
# 1) Read data (results already include outcome_title) ----
# =========================================================
metabo <- readxl::read_xlsx(
  path = "/mnt/project/data/processed/metabolomics_associations.xlsx"
) %>%
  mutate(
    beta_coefficient = as.numeric(gsub(",", ".", beta_coefficient)),
    std_coefficient  = as.numeric(gsub(",", ".", std_coefficient)),
    std_error        = as.numeric(gsub(",", ".", std_error)),
    t_statistic      = as.numeric(gsub(",", ".", t_statistic)),
    p_value_raw      = as.numeric(gsub(",", ".", p_value_raw)),
    p_value_fdr      = as.numeric(gsub(",", ".", p_value_fdr)),
    CI_low           = as.numeric(gsub(",", ".", CI_low)),
    CI_high          = as.numeric(gsub(",", ".", CI_high)),
    significant_fdr  = as.logical(significant_fdr)
  )

# =========================================================
# 2) Create metabolite name 
# =========================================================
metabo_clean <- metabo %>%
  mutate(
    title_use  = dplyr::coalesce(outcome_title, outcome),  # fallback if outcome_title is missing
    metabolite = janitor::make_clean_names(title_use)
  ) %>%
  filter(!is.na(metabolite) & metabolite != "")

# =========================================================
# 3) Define met_classes list
# =========================================================
new_met_list <- c(
  # 1 — CHOLESTEROL (C)
  "Total.C"           = "Total cholesterol",
  "Total.C.noHDL"     = "Total cholesterol minus HDL-C",
  "Remnant.C"         = "Remnant cholesterol (non-HDL, non-LDL -cholesterol)",
  "VLDL.C"            = "VLDL cholesterol",
  "IDL.C"             = "Cholesterol in IDL",
  "LDL.C.clinical"    = "Clinical LDL cholesterol",
  "LDL.C"             = "LDL cholesterol",
  "HDL.C"             = "HDL cholesterol",
  "XXL.VLDL.C"        = "Cholesterol in chylomicrons and extremely large VLDL",
  "XL.VLDL.C"         = "Cholesterol in very large VLDL",
  "L.VLDL.C"          = "Cholesterol in large VLDL",
  "M.VLDL.C"          = "Cholesterol in medium VLDL",
  "S.VLDL.C"          = "Cholesterol in small VLDL",
  "XS.VLDL.C"         = "Cholesterol in very small VLDL",
  "L.LDL.C"           = "Cholesterol in large LDL",
  "M.LDL.C"           = "Cholesterol in medium LDL",
  "S.LDL.C"           = "Cholesterol in small LDL",
  "XL.HDL.C"          = "Cholesterol in very large HDL",
  "L.HDL.C"           = "Cholesterol in large HDL",
  "M.HDL.C"           = "Cholesterol in medium HDL",
  "S.HDL.C"           = "Cholesterol in small HDL",
  "XXL.VLDL.C(%)"     = "Cholesterol to total lipids ratio in chylomicrons and extremely large VLDL",
  "XL.VLDL.C(%)"     = "Cholesterol to total lipids ratio in very large VLDL",
  "L.VLDL.C(%)"       = "Cholesterol to total lipids ratio in large VLDL",
  "M.VLDL.C(%)"       = "Cholesterol to total lipids ratio in medium VLDL",
  "S.VLDL.C(%)"       = "Cholesterol to total lipids ratio in small VLDL",
  "XS.VLDL.C(%)"       = "Cholesterol to total lipids ratio in very small VLDL",
  "IDL.C(%)"       = "Cholesterol to total lipids ratio in IDL",
  "L.LDL.C(%)"     = "Cholesterol to total lipids ratio in large LDL",
  "M.LDL.C(%)"     = "Cholesterol to total lipids ratio in medium LDL",
  "S.LDL.C(%)"     = "Cholesterol to total lipids ratio in small LDL",
  "XL.HDL.C(%)"    = "Cholesterol to total lipids ratio in very large HDL",
  "L.HDL.C(%)"     = "Cholesterol to total lipids ratio in large HDL",
  "M.HDL.C(%)"     = "Cholesterol to total lipids ratio in medium HDL",
  "S.HDL.C(%)"     = "Cholesterol to total lipids ratio in small HDL",
  
  # 2 — TRIGLYCERIDES (TG)
  "Total.TG"        = "Total triglycerides",
  "VLDL.TG"         = "Triglycerides in VLDL",
  "IDL.TG"          = "Triglycerides in IDL",
  "LDL.TG"          = "Triglycerides in LDL",
  "HDL.TG"          = "Triglycerides in HDL",
  
  "XXL.VLDL.TG"     = "Triglycerides in chylomicrons and extremely large VLDL",
  "XL.VLDL.TG"      = "Triglycerides in very large VLDL",
  "L.VLDL.TG"       = "Triglycerides in large VLDL",
  "M.VLDL.TG"       = "Triglycerides in medium VLDL",
  "S.VLDL.TG"       = "Triglycerides in small VLDL",
  "XS.VLDL.TG"      = "Triglycerides in very small VLDL",
  
  "L.LDL.TG"        = "Triglycerides in large LDL",
  "M.LDL.TG"        = "Triglycerides in medium LDL",
  "S.LDL.TG"        = "Triglycerides in small LDL",
  
  "XL.HDL.TG"       = "Triglycerides in very large HDL",
  "L.HDL.TG"        = "Triglycerides in large HDL",
  "M.HDL.TG"        = "Triglycerides in medium HDL",
  "S.HDL.TG"        = "Triglycerides in small HDL",
  
  "XXL.VLDL.TG(%)"  = "Triglycerides to total lipids ratio in chylomicrons and extremely large VLDL",
  "XL.VLDL.TG(%)"   = "Triglycerides to total lipids ratio in very large VLDL",
  "L.VLDL.TG(%)"    = "Triglycerides to total lipids ratio in large VLDL",
  "M.VLDL.TG(%)"    = "Triglycerides to total lipids ratio in medium VLDL",
  "S.VLDL.TG(%)"    = "Triglycerides to total lipids ratio in small VLDL",
  "XS.VLDL.TG(%)"   = "Triglycerides to total lipids ratio in very small VLDL",
  
  "IDL.TG(%)"       = "Triglycerides to total lipids ratio in IDL",
  
  "L.LDL.TG(%)"     = "Triglycerides to total lipids ratio in large LDL",
  "M.LDL.TG(%)"     = "Triglycerides to total lipids ratio in medium LDL",
  "S.LDL.TG(%)"     = "Triglycerides to total lipids ratio in small LDL",
  
  "XL.HDL.TG(%)"    = "Triglycerides to total lipids ratio in very large HDL",
  "L.HDL.TG(%)"     = "Triglycerides to total lipids ratio in large HDL",
  "M.HDL.TG(%)"     = "Triglycerides to total lipids ratio in medium HDL",
  "S.HDL.TG(%)"     = "Triglycerides to total lipids ratio in small HDL",
  
  # 3 — PHOSPHOLIPIDS (PL)
  "Total.PL"        = "Total phospholipids in lipoprotein particles",
  "VLDL.PL"         = "Phospholipids in VLDL",
  "IDL.PL"          = "Phospholipids in IDL",
  "LDL.PL"          = "Phospholipids in LDL",
  "HDL.PL"          = "Phospholipids in HDL",
  "XXL.VLDL.PL"     = "Phospholipids in chylomicrons and extremely large VLDL",
  "XL.VLDL.PL"      = "Phospholipids in very large VLDL",
  "L.VLDL.PL"       = "Phospholipids in large VLDL",
  "M.VLDL.PL"       = "Phospholipids in medium VLDL",
  "S.VLDL.PL"       = "Phospholipids in small VLDL",
  "XS.VLDL.PL"      = "Phospholipids in very small VLDL",
  "L.LDL.PL"        = "Phospholipids in large LDL",
  "M.LDL.PL"        = "Phospholipids in medium LDL",
  "S.LDL.PL"        = "Phospholipids in small LDL",
  "XL.HDL.PL"       = "Phospholipids in very large HDL",
  "L.HDL.PL"        = "Phospholipids in large HDL",
  "M.HDL.PL"        = "Phospholipids in medium HDL",
  "S.HDL.PL"        = "Phospholipids in small HDL",
  "XXL.VLDL.PL(%)"  = "Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL",
  "XL.VLDL.PL(%)"   = "Phospholipids to total lipids ratio in very large VLDL",
  "L.VLDL.PL(%)"    = "Phospholipids to total lipids ratio in large VLDL",
  "M.VLDL.PL(%)"    = "Phospholipids to total lipids ratio in medium VLDL",
  "S.VLDL.PL(%)"    = "Phospholipids to total lipids ratio in small VLDL",
  "XS.VLDL.PL(%)"   = "Phospholipids to total lipids ratio in very small VLDL",
  "IDL.PL(%)"       = "Phospholipids to total lipids ratio in IDL",
  "L.LDL.PL(%)"     = "Phospholipids to total lipids ratio in large LDL",
  "M.LDL.PL(%)"     = "Phospholipids to total lipids ratio in medium LDL",
  "S.LDL.PL(%)"     = "Phospholipids to total lipids ratio in small LDL",
  "XL.HDL.PL(%)"    = "Phospholipids to total lipids ratio in very large HDL",
  "L.HDL.PL(%)"     = "Phospholipids to total lipids ratio in large HDL",
  "M.HDL.PL(%)"     = "Phospholipids to total lipids ratio in medium HDL",
  "S.HDL.PL(%)"     = "Phospholipids to total lipids ratio in small HDL",
  
  # 4 — CHOLESTERYL ESTERS (CE)
  "Total.CE"        = "Total esterified cholesterol",
  "VLDL.CE"         = "Cholesteryl esters in VLDL",
  "IDL.CE"          = "Cholesteryl esters in IDL",
  "LDL.CE"          = "Cholesteryl esters in LDL",
  "HDL.CE"          = "Cholesteryl esters in HDL",
  "XXL.VLDL.CE"     = "Cholesteryl esters in chylomicrons and extremely large VLDL",
  "XL.VLDL.CE"      = "Cholesteryl esters in very large VLDL",
  "L.VLDL.CE"       = "Cholesteryl esters in large VLDL",
  "M.VLDL.CE"       = "Cholesteryl esters in medium VLDL",
  "S.VLDL.CE"       = "Cholesteryl esters in small VLDL",
  "XS.VLDL.CE"      = "Cholesteryl esters in very small VLDL",
  "L.LDL.CE"        = "Cholesteryl esters in large LDL",
  "M.LDL.CE"        = "Cholesteryl esters in medium LDL",
  "S.LDL.CE"        = "Cholesteryl esters in small LDL",
  "XL.HDL.CE"       = "Cholesteryl esters in very large HDL",
  "L.HDL.CE"        = "Cholesteryl esters in large HDL",
  "M.HDL.CE"        = "Cholesteryl esters in medium HDL",
  "S.HDL.CE"        = "Cholesteryl esters in small HDL",
  "XXL.VLDL.CE(%)"  = "Cholesteryl esters to total lipids ratio in chylomicrons and extremely large VLDL",
  "XL.VLDL.CE(%)"   = "Cholesteryl esters to total lipids ratio in very large VLDL",
  "L.VLDL.CE(%)"    = "Cholesteryl esters to total lipids ratio in large VLDL",
  "M.VLDL.CE(%)"    = "Cholesteryl esters to total lipids ratio in medium VLDL",
  "S.VLDL.CE(%)"    = "Cholesteryl esters to total lipids ratio in small VLDL",
  "XS.VLDL.CE(%)"   = "Cholesteryl esters to total lipids ratio in very small VLDL",
  "IDL.CE(%)"       = "Cholesteryl esters to total lipids ratio in IDL",
  "L.LDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in large LDL",
  "M.LDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in medium LDL",
  "S.LDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in small LDL",
  "XL.HDL.CE(%)"    = "Cholesteryl esters to total lipids ratio in very large HDL",
  "L.HDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in large HDL",
  "M.HDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in medium HDL",
  "S.HDL.CE(%)"     = "Cholesteryl esters to total lipids ratio in small HDL",
  
  # 5 — FREE CHOLESTEROL (FC)
  "Total.FC"              = "Total free cholesterol",
  "VLDL.FC"         = "Free cholesterol in VLDL",
  "IDL.FC"          = "Free cholesterol in IDL",
  "LDL.FC"          = "Free cholesterol in LDL",
  "HDL.FC"          = "Free cholesterol in HDL",
  "XXL.VLDL.FC"     = "Free cholesterol in chylomicrons and extremely large VLDL",
  "XL.VLDL.FC"      = "Free cholesterol in very large VLDL",
  "L.VLDL.FC"       = "Free cholesterol in large VLDL",
  "M.VLDL.FC"       = "Free cholesterol in medium VLDL",
  "S.VLDL.FC"       = "Free cholesterol in small VLDL",
  "XS.VLDL.FC"      = "Free cholesterol in very small VLDL",
  "L.LDL.FC"        = "Free cholesterol in large LDL",
  "M.LDL.FC"        = "Free cholesterol in medium LDL",
  "S.LDL.FC"        = "Free cholesterol in small LDL",
  "XL.HDL.FC"       = "Free cholesterol in very large HDL",
  "L.HDL.FC"        = "Free cholesterol in large HDL",
  "M.HDL.FC"        = "Free cholesterol in medium HDL",
  "S.HDL.FC"        = "Free cholesterol in small HDL",
  "XXL.VLDL.FC(%)"  = "Free cholesterol to total lipids ratio in chylomicrons and extremely large VLDL",
  "XL.VLDL.FC(%)"   = "Free cholesterol to total lipids ratio in very large VLDL",
  "L.VLDL.FC(%)"    = "Free cholesterol to total lipids ratio in large VLDL",
  "M.VLDL.FC(%)"    = "Free cholesterol to total lipids ratio in medium VLDL",
  "S.VLDL.FC(%)"    = "Free cholesterol to total lipids ratio in small VLDL",
  "XS.VLDL.FC(%)"   = "Free cholesterol to total lipids ratio in very small VLDL",
  "IDL.FC(%)"       = "Free cholesterol to total lipids ratio in IDL",
  "L.LDL.FC(%)"     = "Free cholesterol to total lipids ratio in large LDL",
  "M.LDL.FC(%)"     = "Free cholesterol to total lipids ratio in medium LDL",
  "S.LDL.FC(%)"     = "Free cholesterol to total lipids ratio in small LDL",
  "XL.HDL.FC(%)"    = "Free cholesterol to total lipids ratio in very large HDL",
  "L.HDL.FC(%)"     = "Free cholesterol to total lipids ratio in large HDL",
  "M.HDL.FC(%)"     = "Free cholesterol to total lipids ratio in medium HDL",
  "S.HDL.FC(%)"     = "Free cholesterol to total lipids ratio in small HDL",
  
  # 6 — TOTAL LIPIDS (TL)
  "Total lipids"        = "Total lipids in lipoprotein particles",
  "VLDL.Total lipids"   = "Total lipids in VLDL",
  "IDL.Total lipids"    = "Total lipids in IDL",
  "LDL.Total lipids"    = "Total lipids in LDL",
  "HDL.Total lipids"    = "Total lipids in HDL",
  "XXL.VLDL.Total lipids"   = "Total lipids in chylomicrons and extremely large VLDL",
  "XL.VLDL.Total lipids"    = "Total lipids in very large VLDL",
  "L.VLDL.Total lipids"     = "Total lipids in large VLDL",
  "M.VLDL.Total lipids"     = "Total lipids in medium VLDL",
  "S.VLDL.Total lipids"     = "Total lipids in small VLDL",
  "XS.VLDL.Total lipids"    = "Total lipids in very small VLDL",
  "L.LDL.Total lipids"      = "Total lipids in large LDL",
  "M.LDL.Total lipids"      = "Total lipids in medium LDL",
  "S.LDL.Total lipids"      = "Total lipids in small LDL",
  "XL.HDL.Total lipids"     = "Total lipids in very large HDL",
  "L.HDL.Total lipids"       = "Total lipids in large HDL",
  "M.HDL.Total lipids"      = "Total lipids in medium HDL",
  "S.HDL.Total lipids"      = "Total lipids in small HDL",
  
  # 7 — LIPOPROTEIN PARTICLE CONCENTRATIONS (LPC)
  "Total.lipoprotein.concentration"        = "Total concentration of lipoprotein particles",
  "VLDL.concentration"         = "Concentration of VLDL particles",
  "IDL.concentration"          = "Concentration of IDL particles",
  "LDL.concentration"          = "Concentration of LDL particles",
  "HDL.concentration"          = "Concentration of HDL particles",
  "XXL.VLDL.concentration"     = "Concentration of chylomicrons and extremely large VLDL particles",
  "XL.VLDL.concentration"      = "Concentration of very large VLDL particles",
  "L.VLDL.concentration"       = "Concentration of large VLDL particles",
  "M.VLDL.concentration"       = "Concentration of medium VLDL particles",
  "S.VLDL.concentration"       = "Concentration of small VLDL particles",
  "XS.VLDL.concentration"      = "Concentration of very small VLDL particles",
  "L.LDL.concentration"        = "Concentration of large LDL particles",
  "M.LDL.concentration"        = "Concentration of medium LDL particles",
  "S.LDL.concentration"        = "Concentration of small LDL particles",
  "XL.HDL.concentration"       = "Concentration of very large HDL particles",
  "L.HDL.concentration"        = "Concentration of large HDL particles",
  "M.HDL.concentration"        = "Concentration of medium HDL particles",
  "S.HDL.concentration"        = "Concentration of small HDL particles",
  
  # 8 — PHOSPHOLIPID-RELATED METABOLITES
  "Total cholines"       = "Total cholines",
  "Phosphatidylcholines" = "Phosphatidylcholines",
  "Phosphoglycerides"    = "Phosphoglycerides",
  "TG.Phosphoglycerides(%)" = "Ratio of triglycerides to phosphoglycerides",
  "Sphingomyelins"          = "Sphingomyelins",
  
  # 9 — LIPOPROTEIN PARTICLE SIZES (LPS)
  "VLDL.size"  = "Average diameter for VLDL particles",
  "LDL.size"   = "Average diameter for LDL particles",
  "HDL.size"   = "Average diameter for HDL particles",
  
  # 10 — FATTY ACIDS (FA)
  "Total.FA"          = "Total fatty acids",
  "Unsaturation.FA"   = "Degree of unsaturation",
  "Omega-3"         = "Omega-3 fatty acids",
  "Omega-6"         = "Omega-6 fatty acids",
  "Pufa.FA"           = "Polyunsaturated fatty acids",
  "Mufa.FA"           = "Monounsaturated fatty acids",
  "Saturated.FA"            = "Saturated fatty acids",
  "Linoleic.Acid"             = "Linoleic acid",
  "Docosahexaenoic.Acid"            = "Docosahexaenoic acid",
  "Omega-3.FA(%)"      = "Ratio of omega-3 fatty acids to total fatty acids",
  "Omega-6.FA(%)"      = "Ratio of omega-6 fatty acids to total fatty acids",
  "Pufa.FA(%)"         = "Ratio of polyunsaturated fatty acids to total fatty acids",
  "Mufa.FA(%)"        = "Ratio of monounsaturated fatty acids to total fatty acids",
  "Saturated.FA(%)"         = "Ratio of saturated fatty acids to total fatty acids",
  "Linoleic.Acid.FA(%)"          = "Ratio of linoleic acid to total fatty acids",
  "Docosahexaenoic.Acid.FA(%)"         = "Ratio of docosahexaenoic acid to total fatty acids",
  "Pufa.Mufa.FA(%)" = "Ratio of polyunsaturated fatty acids to monounsaturated fatty acids",
  "Omega-6.Omega-3(%)"  = "Ratio of omega-6 fatty acids to omega-3 fatty acids",
  
  # 11 — AMINO ACIDS (AA)
  "Alanine"        = "Alanine",
  "Glutamine"      = "Glutamine",
  "Glycine"        = "Glycine",
  "Histidine"      = "Histidine",
  "BCAA.concentration"           = "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)",
  "Isoleucine"     = "Isoleucine",
  "Leucine"        = "Leucine",
  "Valine"         = "Valine",
  "Phenylalanine"  = "Phenylalanine",
  "Tyrosine"       = "Tyrosine",
  "Alanine.Corrected" = "Spectrometer-corrected alanine",
  
  # 12 — GLYCOLYSIS (GLY)
  "Glucose"   = "Glucose",
  "Lactate"   = "Lactate",
  "Pyruvate"       = "Pyruvate",
  "Citrate"        = "Citrate",
  "Glucose.Lactate"        = "Glucose-lactate",
  
  # 13 — KETONE BODIES (KET)
  "3-Hydroxybutyrate"    = "3-Hydroxybutyrate",
  "Acetate"       = "Acetate",
  "Acetoacetate"  = "Acetoacetate",
  "Acetone"       = "Acetone",
  
  # 14 — FLUID BALANCE (FB)
  "Creatinine"   = "Creatinine",
  "Albumin"      = "Albumin",
  
  # 15 — Inflammation
  "GlycA" = "Glycoprotein acetyls",
  
  # 16 — APOLIPOPROTEINS (APO)
  "APO_B"        = "Apolipoprotein B",
  "APO_A1"       = "Apolipoprotein A1",
  "APO_B.APO_A1(%)" = "Ratio of apolipoprotein B to apolipoprotein A1"
)

length(new_met_list)  # it should be 251 metabolites


stopifnot(exists("new_met_list"))
stopifnot(length(new_met_list) == 251)

# =========================================================
# 3b) Define 16 functional families (short-name based) 
# =========================================================
family_map <- list(
  Cholesterol = c(
    "Total.C","Total.C.noHDL","Remnant.C","VLDL.C","IDL.C","LDL.C.clinical","LDL.C","HDL.C",
    "XXL.VLDL.C","XL.VLDL.C","L.VLDL.C","M.VLDL.C","S.VLDL.C","XS.VLDL.C",
    "L.LDL.C","M.LDL.C","S.LDL.C",
    "XL.HDL.C","L.HDL.C","M.HDL.C","S.HDL.C",
    "XXL.VLDL.C(%)","XL.VLDL.C(%)","L.VLDL.C(%)","M.VLDL.C(%)","S.VLDL.C(%)","XS.VLDL.C(%)",
    "IDL.C(%)",
    "L.LDL.C(%)","M.LDL.C(%)","S.LDL.C(%)",
    "XL.HDL.C(%)","L.HDL.C(%)","M.HDL.C(%)","S.HDL.C(%)"
  ),
  Triglycerides = c(
    "Total.TG","VLDL.TG","IDL.TG","LDL.TG","HDL.TG",
    "XXL.VLDL.TG","XL.VLDL.TG","L.VLDL.TG","M.VLDL.TG","S.VLDL.TG","XS.VLDL.TG",
    "L.LDL.TG","M.LDL.TG","S.LDL.TG",
    "XL.HDL.TG","L.HDL.TG","M.HDL.TG","S.HDL.TG",
    "XXL.VLDL.TG(%)","XL.VLDL.TG(%)","L.VLDL.TG(%)","M.VLDL.TG(%)","S.VLDL.TG(%)","XS.VLDL.TG(%)",
    "IDL.TG(%)",
    "L.LDL.TG(%)","M.LDL.TG(%)","S.LDL.TG(%)",
    "XL.HDL.TG(%)","L.HDL.TG(%)","M.HDL.TG(%)","S.HDL.TG(%)"
  ),
  Phospholipids = c(
    "Total.PL","VLDL.PL","IDL.PL","LDL.PL","HDL.PL",
    "XXL.VLDL.PL","XL.VLDL.PL","L.VLDL.PL","M.VLDL.PL","S.VLDL.PL","XS.VLDL.PL",
    "L.LDL.PL","M.LDL.PL","S.LDL.PL",
    "XL.HDL.PL","L.HDL.PL","M.HDL.PL","S.HDL.PL",
    "XXL.VLDL.PL(%)","XL.VLDL.PL(%)","L.VLDL.PL(%)","M.VLDL.PL(%)","S.VLDL.PL(%)","XS.VLDL.PL(%)",
    "IDL.PL(%)",
    "L.LDL.PL(%)","M.LDL.PL(%)","S.LDL.PL(%)",
    "XL.HDL.PL(%)","L.HDL.PL(%)","M.HDL.PL(%)","S.HDL.PL(%)"
  ),
  Cholesteryl_esters = c(
    "Total.CE","VLDL.CE","IDL.CE","LDL.CE","HDL.CE",
    "XXL.VLDL.CE","XL.VLDL.CE","L.VLDL.CE","M.VLDL.CE","S.VLDL.CE","XS.VLDL.CE",
    "L.LDL.CE","M.LDL.CE","S.LDL.CE",
    "XL.HDL.CE","L.HDL.CE","M.HDL.CE","S.HDL.CE",
    "XXL.VLDL.CE(%)","XL.VLDL.CE(%)","L.VLDL.CE(%)","M.VLDL.CE(%)","S.VLDL.CE(%)","XS.VLDL.CE(%)",
    "IDL.CE(%)",
    "L.LDL.CE(%)","M.LDL.CE(%)","S.LDL.CE(%)",
    "XL.HDL.CE(%)","L.HDL.CE(%)","M.HDL.CE(%)","S.HDL.CE(%)"
  ),
  Free_cholesterol = c(
    "Total.FC","VLDL.FC","IDL.FC","LDL.FC","HDL.FC",
    "XXL.VLDL.FC","XL.VLDL.FC","L.VLDL.FC","M.VLDL.FC","S.VLDL.FC","XS.VLDL.FC",
    "L.LDL.FC","M.LDL.FC","S.LDL.FC",
    "XL.HDL.FC","L.HDL.FC","M.HDL.FC","S.HDL.FC",
    "XXL.VLDL.FC(%)","XL.VLDL.FC(%)","L.VLDL.FC(%)","M.VLDL.FC(%)","S.VLDL.FC(%)","XS.VLDL.FC(%)",
    "IDL.FC(%)",
    "L.LDL.FC(%)","M.LDL.FC(%)","S.LDL.FC(%)",
    "XL.HDL.FC(%)","L.HDL.FC(%)","M.HDL.FC(%)","S.HDL.FC(%)"
  ),
  Total_lipids = c(
    "Total lipids","VLDL.Total lipids","IDL.Total lipids","LDL.Total lipids","HDL.Total lipids",
    "XXL.VLDL.Total lipids","XL.VLDL.Total lipids","L.VLDL.Total lipids","M.VLDL.Total lipids",
    "S.VLDL.Total lipids","XS.VLDL.Total lipids",
    "L.LDL.Total lipids","M.LDL.Total lipids","S.LDL.Total lipids",
    "XL.HDL.Total lipids","L.HDL.Total lipids","M.HDL.Total lipids","S.HDL.Total lipids"
  ),
  Lipoprotein_particle_conc = c(
    "Total.lipoprotein.concentration","VLDL.concentration","IDL.concentration","LDL.concentration","HDL.concentration",
    "XXL.VLDL.concentration","XL.VLDL.concentration","L.VLDL.concentration","M.VLDL.concentration",
    "S.VLDL.concentration","XS.VLDL.concentration",
    "L.LDL.concentration","M.LDL.concentration","S.LDL.concentration",
    "XL.HDL.concentration","L.HDL.concentration","M.HDL.concentration","S.HDL.concentration"
  ),
  Phospholipid_related = c(
    "Total cholines","Phosphatidylcholines","Phosphoglycerides","TG.Phosphoglycerides(%)","Sphingomyelins"
  ),
  Lipoprotein_particle_sizes = c(
    "VLDL.size","LDL.size","HDL.size"
  ),
  Fatty_acids = c(
    "Total.FA","Unsaturation.FA","Omega-3","Omega-6","Pufa.FA","Mufa.FA","Saturated.FA",
    "Linoleic.Acid","Docosahexaenoic.Acid",
    "Omega-3.FA(%)","Omega-6.FA(%)","Pufa.FA(%)","Mufa.FA(%)","Saturated.FA(%)",
    "Linoleic.Acid.FA(%)","Docosahexaenoic.Acid.FA(%)",
    "Pufa.Mufa.FA(%)","Omega-6.Omega-3(%)"
  ),
  Amino_acids = c(
    "Alanine","Glutamine","Glycine","Histidine","BCAA.concentration",
    "Isoleucine","Leucine","Valine","Phenylalanine","Tyrosine","Alanine.Corrected"
  ),
  Glycolysis = c("Glucose","Lactate","Pyruvate","Citrate","Glucose.Lactate"),
  Ketone_bodies = c("3-Hydroxybutyrate","Acetate","Acetoacetate","Acetone"),
  Fluid_balance = c("Creatinine","Albumin"),
  Inflammation = c("GlycA"),
  Apolipoproteins = c("APO_B","APO_A1","APO_B.APO_A1(%)")
)

# =========================================================
# 4) Build lookup (snake_case title -> short name) and recode
# =========================================================
stopifnot(exists("new_met_list"))
stopifnot(length(new_met_list) == 251)

lookup <- setNames(
  names(new_met_list),                              # short name (e.g., "Total.C")
  janitor::make_clean_names(unname(new_met_list))   # snake_case of official title
)

metabo_final <- metabo_clean %>%
  mutate(Metabolite = unname(lookup[metabolite])) %>%
  filter(!is.na(Metabolite))

# =========================================================
# 5) Attach functional family (short-name based)
# =========================================================
stopifnot(exists("family_map"))

family_df <- do.call(rbind, lapply(names(family_map), function(fam) {
  data.frame(
    Metabolite = family_map[[fam]],
    family = fam,
    stringsAsFactors = FALSE
  )
}))

# (Optional sanity check) Which short-names are not classified?
all_short <- names(new_met_list)
missing_family <- setdiff(all_short, family_df$Metabolite)
if (length(missing_family) > 0) {
  message("⚠️ Unclassified variables (not in family_map): ", paste(missing_family, collapse = ", "))
}

metabo_final <- metabo_final %>%
  left_join(family_df, by = "Metabolite") %>%
  mutate(
    family = if_else(is.na(family), "Unclassified", family),
    ChemicalClass = factor(family, levels = c(names(family_map), "Unclassified"))
  )

# =========================================================
# 6) Keep compatibility names for p-values
# =========================================================
metabo_final <- metabo_final %>%
  rename(
    p.value = p_value_raw,
    p_fdr   = p_value_fdr
  )

# =========================================================
# 7) Filter significant metabolites (FDR) for plotting 
# =========================================================
metabo_sig <- metabo_final %>%
  filter(p_fdr < 0.05,
         !is.na(beta_coefficient),
         abs(beta_coefficient) > 0.01,
         !is.na(Metabolite))

# =========================================================
# 8) Helpers: lipoprotein class/size for structured ordering
# =========================================================
get_lipoprotein_class <- function(x) {
  dplyr::case_when(
    grepl("VLDL", x) ~ "VLDL",
    grepl("IDL",  x) ~ "IDL",
    grepl("LDL",  x) ~ "LDL",
    grepl("HDL",  x) ~ "HDL",
    TRUE ~ "Other"
  )
}

get_lipoprotein_size <- function(x) {
  dplyr::case_when(
    grepl("^XXL\\.", x) ~ "XXL",
    grepl("^XL\\.",  x) ~ "XL",
    grepl("^L\\.",   x) ~ "L",
    grepl("^M\\.",   x) ~ "M",
    grepl("^S\\.",   x) ~ "S",
    grepl("^XS\\.",  x) ~ "XS",
    TRUE ~ "Other"
  )
}

lipoprotein_class_levels <- c("VLDL","IDL","LDL","HDL","Other")
lipoprotein_size_levels  <- c("XXL","XL","L","M","S","XS","Other")

# =========================================================
# 9) Prepare ordered dataframe + spacer rows
# =========================================================
if (nrow(metabo_sig) == 0) stop("No significant metabolites at FDR < 0.05")

metabo_plot_df <- metabo_sig %>%
  mutate(
    LipoproteinClass = factor(get_lipoprotein_class(Metabolite), levels = lipoprotein_class_levels),
    LipoproteinSize  = factor(get_lipoprotein_size(Metabolite),  levels = lipoprotein_size_levels)
  ) %>%
  arrange(ChemicalClass, LipoproteinClass, LipoproteinSize, Metabolite) %>%
  mutate(
    idx = row_number(),
    boundary = !is.na(lag(LipoproteinClass)) &
      (ChemicalClass == lag(ChemicalClass)) &
      (LipoproteinClass != lag(LipoproteinClass))
  )

# Insert invisible bars before each change of LipoproteinClass (within each family)
spacers <- metabo_plot_df %>%
  filter(boundary) %>%
  transmute(
    Metabolite = NA_character_,
    family = family,
    ChemicalClass = ChemicalClass,
    LipoproteinClass = LipoproteinClass,
    LipoproteinSize = factor("Other", levels = lipoprotein_size_levels),
    beta_coefficient = 0,
    Spacer = TRUE,
    pos = idx - 0.5
  )

metabo_plot_df <- metabo_plot_df %>%
  mutate(Spacer = FALSE, pos = idx) %>%
  bind_rows(spacers) %>%
  arrange(pos) %>%
  mutate(
    id = row_number()
  )

# Label offset proportional to the beta range (prevents the “all labels glued” effect)
beta_rng <- range(metabo_plot_df$beta_coefficient, na.rm = TRUE)
delta <- max(0.01, 0.03 * max(abs(beta_rng)))  # adaptive small offset

metabo_plot_df <- metabo_plot_df %>%
  mutate(
    label_y = ifelse(!Spacer,
                     ifelse(beta_coefficient >= 0,
                            beta_coefficient + delta,   
                            pmax(beta_coefficient + delta, 0.01)
                     ),
                     NA_real_)
  )

# =========================================================
# 10) Color palette
# =========================================================
class_cols <- setNames(
  scales::hue_pal()(length(levels(metabo_plot_df$ChemicalClass))),
  levels(metabo_plot_df$ChemicalClass)
)

# =========================================================
# 11) Circular plot
# =========================================================

pad <- 0.005
bmin <- min(metabo_plot_df$beta_coefficient[!metabo_plot_df$Spacer], na.rm = TRUE)
bmax <- max(metabo_plot_df$beta_coefficient[!metabo_plot_df$Spacer], na.rm = TRUE)

# Limits
lim_low  <- floor((bmin - pad) / 0.05) * 0.05
lim_high <- ceiling((bmax + pad) / 0.05) * 0.05

p_gg <- ggplot(metabo_plot_df, aes(x = factor(id), y = beta_coefficient, fill = ChemicalClass)) +
  geom_col(aes(alpha = ifelse(Spacer, 0, 1)), width = 0.98) +
  scale_alpha_identity(guide = "none") +
  ggplot2::coord_radial(
    start = 0,
    end = -0.04 * pi,
    expand = TRUE,
    inner.radius = 0.45,
    rotate.angle = TRUE
  ) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_text(
    data = metabo_plot_df %>% dplyr::filter(!Spacer),
    aes(y = label_y, label = Metabolite, angle = 90, hjust = 0),
    size = 2
  ) +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(lim_low, lim_high),
    breaks = seq(lim_low, lim_high, by = 0.05),
    expand = expansion(mult = c(0, 0))  # evita margen extra raro
  ) +
  labs(title = "Associations between VO2max and metabolites", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(10, 20, 10, 10)
  )

print(p_gg)


# =========================================================
# 12) Export vector PDF ----
# =========================================================
if (!requireNamespace("Cairo", quietly = TRUE)) install.packages("Cairo")

ggsave(
  filename = "./upload/figure_metabolomics.pdf",
  plot = p_gg,
  width = 10,
  height = 8,
  device = "pdf"
)

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
