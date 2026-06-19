# Load necessary libraries
library(readr)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ReactomePA)

### OBSOLETE: This script performs per-miRNA pathway enirhcment. However, many of the individual miRNAs don't have enriched pathways therefore it needs to be performed together.

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold
mirna_count <- 1 # Minimum number of miRNAs required to target a target gene
database_count <- 1 # Minimum number of databases a target gene must appear in


#==========================#
# Create Background Universe 
#==========================#

# Read in signficant results from dose 24 vs control
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Create list of miRNA IDs
dose24_mirnas_all <- dose24_data$Symbol

# Identify gene targets for all DE miRNAs
target_results_background <- get_multimir(mirna = dose24_mirnas_all, table = "mirtarbase")

# Remove redundant target genes
target_results_background_unique <- unique(target_results_background@data$target_entrez) # Entrez format genes


#==========================#
# Identify Target Genes
#==========================#

# Extract significant UP- and DOWN-regulated miRNAs
dose24_data_sig_up <- dose24_data %>% filter(logFC > logFC_threshold) # Upregulated
dose24_data_sig_down <- dose24_data %>% filter(logFC < -logFC_threshold) # Downregulated

# Create list of miRNA IDs
dose24_miRNA_list_up <- dose24_data_sig_up$Symbol
dose24_miRNA_list_down <- dose24_data_sig_down$Symbol

# Identify gene targets for each miRNA
target_results_up <- get_multimir(mirna = dose24_miRNA_list_up, table = "mirtarbase")
target_results_down <- get_multimir(mirna = dose24_miRNA_list_down, table = "mirtarbase")

# Remove redundant target genes
target_results_up_unique <- unique(target_results_up@data$target_entrez)
target_results_down_unique <- unique(target_results_down@data$target_entrez)

#==========================#
# Pathway Enrichment
#==========================#


for (mirna_symbol in dose24_miRNA_list_up) {
  
  targets <- get_multimir(mirna = mirna_symbol, table = 'mirtarbase')
  
  genes <- targets@data$target_entrez
  
  genes_unique <- unique((targets@data$target_entrez))
  
  print(mirna_symbol)
  
  print('Genes:')
  
  print(length(genes_unique))
  
  epath <- enrichPathway(genes_unique, universe = target_results_background_unique)
  
  print('Pathways:')
  
  print(length(epath@result$p.adjust[epath@result$p.adjust < 1]))
  
  print('Pathways below Padj 0.05:')
  
  print(length(epath@result$p.adjust[epath@result$p.adjust < 0.05]))

}







