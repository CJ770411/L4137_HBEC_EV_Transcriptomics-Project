# Load necessary libraries
library(readr)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ReactomePA)

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold for ORA
FDR_threshold <- 0.05 # FDR significance threshold for ORA
set.seed(123) # Set random seed

#==========================#
# Identify Target Genes
#==========================#

# Read in signficant results from dose 24 vs control
dose24_data <- read_csv("../../../results/primaryHBEC/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Extract significant UP- and DOWN-regulated miRNAs
dose24_up <- dose24_data %>% filter(logFC >= logFC_threshold, FDR < FDR_threshold) %>% pull(Symbol) # Upregulated
dose24_down <- dose24_data %>% filter(logFC <= -logFC_threshold, FDR < FDR_threshold) %>% pull(Symbol) # Downregulated

# Identify gene targets for each miRNA
target_results_up <- get_multimir(mirna = dose24_up, table = "mirtarbase")
target_results_down <- get_multimir(mirna = dose24_down, table = "mirtarbase")

# Remove redundant target genes
target_results_up_unique <- unique(target_results_up@data$target_entrez)
target_results_down_unique <- unique(target_results_down@data$target_entrez)

# Print number of target genes identified
length(target_results_up_unique)
length(target_results_down_unique)

#==========================#
# Pathway Enrichment
#==========================#

# Perform pathway enrichment
epath_up <- enrichPathway(target_results_up_unique)
epath_down <- enrichPathway(target_results_down_unique)

# Visualise pathway enrichment
dotplot(epath_up)
dotplot(epath_down)

# View results as a table
epath_up_results_table <- epath_up@result 
epath_down_results_table <- epath_down@result 

# Save enrichment results table as object 
saveRDS(
  epath_up_results_table,
  "../../../results/primaryHBEC/08_FE_indirect/epath_up_results_table.rds"
)

saveRDS(
  epath_down_results_table,
  "../../../results/primaryHBEC/08_FE_indirect/epath_down_results_table.rds"
)

# Save enrichment object to create downstream results figures
saveRDS(
  epath_up,
  "../../../results/primaryHBEC/08_FE_indirect/epath_up_results.rds"
)

saveRDS(
  epath_down,
  "../../../results/primaryHBEC/08_FE_indirect/epath_down_results.rds"
)








