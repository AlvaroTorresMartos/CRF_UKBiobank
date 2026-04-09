
# Useful resources -----
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# https://rpubs.com/pranali018/enrichmentAnalysis

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(clusterProfiler, enrichplot, dplyr, org.Hs.eg.db, ggplot2, stringr)

# KEGG only significant ----
## Up-regulated -----
proteomics = readxl::read_xlsx("./proteomics_associations.xlsx")

proteomics = proteomics %>%
  dplyr::mutate(outcome = stringr::str_to_upper(outcome)) 

universe = proteomics$outcome

proteomics = proteomics %>% 
  dplyr::filter(p_value_fdr < 0.05) %>% 
  dplyr::filter(beta_coefficient > 0)

keytypes(org.Hs.eg.db)


converted = bitr(proteomics$outcome, 
                 fromType = "ALIAS", 
                 toType = "ENTREZID", 
                 OrgDb = 'org.Hs.eg.db')

universe = bitr(universe,
                fromType = "ALIAS",
                toType = "ENTREZID",
                OrgDb = 'org.Hs.eg.db')

ego_kegg = enrichKEGG(gene  = converted$ENTREZID,
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 # universe = universe$ENTREZID,
                 minGSSize = 1, 
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr"
                 )

head(ego_kegg, n = 20)

ego_kegg = readRDS("./ORA_kegg_downregulated.RDS")

x = setReadable(ego_kegg, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)

results = x2@result %>% 
  dplyr::arrange(p.adjust) %>% 
  head(n = 30) 

results$Description %>% head(n = 30)

# saveRDS(ego_kegg, "./ORA_kegg_upregulated.RDS")

# writexl::write_xlsx(results, "./ora_kegg_upregulated.xlsx")


dotplot(x2, showCategory = c(
  # "Cell adhesion molecules",
  # "PI3K-Akt signaling pathway",
  "Cytoskeleton in muscle cells",
  # "Cytokine-cytokine receptor interaction",
  "ECM-receptor interaction",
  "Pancreatic secretion",
                                   "Nitrogen metabolism", 
                             "Focal adhesion", 
                             "Ras signaling pathway", 
                             "Hormone signaling",                              
                             "Fat digestion and absorption", 
                             "Regulation of actin cytoskeleton", 
                             # "MAPK signaling pathway",                          
                             "Cholesterol metabolism", 
                             "PPAR signaling pathway"                           
                             # "JAK-STAT signaling pathway"   
                                   ))


# upsetplot(x2)
# heatplot(x2)
cnetplot(x2, showCategory = c(# "Cell adhesion molecules",
  # "PI3K-Akt signaling pathway",
  "Cytoskeleton in muscle cells",
  # "Cytokine-cytokine receptor interaction",
  "ECM-receptor interaction",
  "Pancreatic secretion",
  "Nitrogen metabolism",
  "Focal adhesion",
  "Ras signaling pathway", 
  "Hormone signaling",                              
  "Fat digestion and absorption", 
  "Regulation of actin cytoskeleton", 
  # "MAPK signaling pathway",                          
  "Cholesterol metabolism", 
  "PPAR signaling pathway"                        
  # "JAK-STAT signaling pathway" 
  ))
# emapplot(x2)
# treeplot(x2)

## Down-regulated -----
proteomics = readxl::read_xlsx("./proteomics_associations.xlsx")

proteomics = proteomics %>%
  dplyr::mutate(outcome = stringr::str_to_upper(outcome)) 

universe = proteomics$outcome

proteomics = proteomics %>% 
  dplyr::filter(p_value_fdr < 0.05) %>% 
  dplyr::filter(beta_coefficient < 0)

keytypes(org.Hs.eg.db)


converted = bitr(proteomics$outcome, 
                 fromType = "ALIAS", 
                 toType = "ENTREZID", 
                 OrgDb = 'org.Hs.eg.db')

universe = bitr(universe,
                fromType = "ALIAS",
                toType = "ENTREZID",
                OrgDb = 'org.Hs.eg.db')

ego_kegg = enrichKEGG(gene  = converted$ENTREZID,
                      keyType = 'ncbi-geneid',
                      organism = 'hsa',
                      # universe = universe$ENTREZID,
                      minGSSize = 1, 
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr"
)

ego_kegg = readRDS("./ORA_kegg_downregulated.RDS")

x = setReadable(ego_kegg, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)

results = x2@result %>% 
  dplyr::arrange(p.adjust) %>% 
  head (n = 30)

results$Description %>% head(n = 30)

# saveRDS(ego_kegg, "./ORA_kegg_downregulated.RDS")

writexl::write_xlsx(results, "./ora_kegg_downregulated.xlsx")

dotplot(x2, showCategory = c(
  "Cytokine-cytokine receptor interaction",
  "Complement and coagulation cascades", 
  "Lipid and atherosclerosis",                                   
  "Apoptosis", 
  "TNF signaling pathway", 
  "NF-kappa B signaling pathway",                                 
  "B cell receptor signaling pathway",
  "MAPK signaling pathway",                                     
  "Cell adhesion molecules",                                      
  "PI3K-Akt signaling pathway",                                   
  "IL-17 signaling pathway" ,                                    
  "Chemokine signaling pathway",                                 
  "AGE-RAGE signaling pathway in diabetic complications",       
  "Toll-like receptor signaling pathway",                        
  "JAK-STAT signaling pathway",                                   
  "HIF-1 signaling pathway",                                    
  "Fluid shear stress and atherosclerosis"
))


# upsetplot(x2)
# heatplot(x2)
cnetplot(x2, showCategory = c("Cytokine-cytokine receptor interaction",
                              "PI3K-Akt signaling pathway",
                              "MAPK signaling pathway"#,
                              # "Lipid and atherosclerosis"
                              # "Complement and coagulation cascades", 
                              #                                   
                              # "Apoptosis", 
                              # "TNF signaling pathway", 
                              # "NF-kappa B signaling pathway",                                 
                              # "B cell receptor signaling pathway",
                             #,                                     
                              # "Cell adhesion molecules",                                      
                              #                                    
                              # "IL-17 signaling pathway" ,                                    
                              # "Chemokine signaling pathway",                                 
                              # "AGE-RAGE signaling pathway in diabetic complications",       
                              # "Toll-like receptor signaling pathway",                        
                              # "JAK-STAT signaling pathway",                                   
                              # "HIF-1 signaling pathway",                                    
                              # "Fluid shear stress and atherosclerosis"
                              ))
# emapplot(x2)
# treeplot(x2)

