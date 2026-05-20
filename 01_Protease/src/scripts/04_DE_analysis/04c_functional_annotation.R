#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("multiMiR")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")

# Load necessary libraries
library(readr)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)

dose24_data <- read_csv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv")

dose24_miRNA_list <- dose24_data$Symbol

target_results <- get_multimir(mirna = dose24_miRNA_list, table = "validated")

genes <- unique(target_results@data$target_symbol)


ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
dotplot(ego)

genes_entrezid <- bitr(genes, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ekegg <- enrichKEGG(genes_entrezid$ENTREZID)

dotplot(ekegg)

