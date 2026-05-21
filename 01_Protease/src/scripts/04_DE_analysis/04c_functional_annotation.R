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

# Read in results from '04b_DE_analysis_edgeR'
dose24_data <- read_csv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv")

# Create list of miRNA IDs
dose24_miRNA_list <- dose24_data$Symbol

# Identify gene targets for each miRNA
target_results <- get_multimir(mirna = dose24_miRNA_list, table = "validated")

# Extract target gene symbols (for GO enrichment compatibility)
genes_symbol <- unique(target_results@data$target_symbol)

# Perform GO enrichment analysis 
ego <- enrichGO(gene = genes_symbol,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# Visualise GO enrichment analysis results
dotplot(ego)

# Extract target gene entrez ID (for KEGG enrichment compatibility)
genes_entrez <- unique(target_results@data$target_entrez)

# Perform KEGG enrichment analysis
ekegg <- enrichKEGG(genes_entrez)

# Visualise KEGG enrichment analysis results
dotplot(ekegg)

# [Outstanding: miRNA interaction network w/ Cytoscope]

