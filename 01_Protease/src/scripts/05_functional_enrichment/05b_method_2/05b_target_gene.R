# Load necessary libraries
library(readr)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)


# Read in results from '04b_DE_analysis_edgeR'
dose24_data <- read_csv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv")

# Extract UP- and DOWN-regulated miRNAs
dose24_data_up <- dose24_data %>% filter(logFC > 0) # Upregulated
dose24_data_down <- dose24_data %>% filter(logFC < 0) # Downregulated

# Create list of miRNA IDs
dose24_miRNA_list_up <- dose24_data_up$Symbol
dose24_miRNA_list_down <- dose24_data_down$Symbol

# Identify gene targets for each miRNA
target_results_up <- get_multimir(mirna = dose24_miRNA_list_up, table = "validated")
target_results_down <- get_multimir(mirna = dose24_miRNA_list_down, table = "validated")

# Extract target gene symbols (for GO enrichment compatibility)
genes_symbol_up <- unique(target_results_up@data$target_symbol)
genes_symbol_down <- unique(target_results_down@data$target_symbol)

# Perform GO enrichment analysis 
ego_up <- enrichGO(gene = genes_symbol_up,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

ego_down <- enrichGO(gene = genes_symbol_down,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

# Visualise GO enrichment analysis results
dotplot(ego_up)
dotplot(ego_down)

# Extract target gene entrez ID (for KEGG enrichment compatibility)
genes_entrez_up <- unique(target_results_up@data$target_entrez)
genes_entrez_down <- unique(target_results_down@data$target_entrez)

# Perform KEGG enrichment analysis
ekegg_up <- enrichKEGG(genes_entrez_up)
ekegg_down <- enrichKEGG(genes_entrez_down)

# Visualise KEGG enrichment analysis results
dotplot(ekegg_up)
dotplot(ekegg_down)

# Test code (not part of main script)
target_results_up@data %>% filter(support_type != "Functional MTI (Weak)")

target_results_up@data %>%
  count(support_type)




# [split gene list by upregulated and downregulated mirnas]
# [Outstanding: 
#              -  enrichPathway()]
#              -  miRNA interaction network w/ Cytoscope]

