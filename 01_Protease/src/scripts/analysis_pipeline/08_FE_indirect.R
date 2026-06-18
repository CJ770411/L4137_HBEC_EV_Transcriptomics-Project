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

logFC_threshold <- 1 # miRNA expression threshold
mirna_count <- 1 # Minimum number of miRNAs required to target a target gene
database_count <- 1 # Minimum number of databases a target gene must appear in


#==========================#
# Create Background Universe 
#==========================#

# [] needs to have separate background universes for upregulated and downregulated de mirnas

# Read in signficant results from dose 24 vs control
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Extract UP- and DOWN-regulated miRNAs
dose24_mirnas_all_up <- dose24_data %>% filter(logFC > 0) # Upregulated
dose24_mirnas_all_down <- dose24_data %>% filter(logFC < 0) # Downregulated

# Create list of miRNA IDs
mirnas_all_up <- dose24_mirnas_all_up$Symbol
mirnas_all_down <- dose24_mirnas_all_down$Symbol

# Identify gene targets for all DE miRNAs
target_results_background_up <- get_multimir(mirna = mirnas_all_up, table = "validated")
target_results_background_down <- get_multimir(mirna = mirnas_all_down, table = "validated")

# Remove redundant target genes
target_results_background_up_unique <- unique(target_results_background_up@data$target_entrez)
target_results_background_down_unique <- unique(target_results_background_down@data$target_entrez)


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
target_results_up <- get_multimir(mirna = dose24_miRNA_list_up, table = "validated")
target_results_down <- get_multimir(mirna = dose24_miRNA_list_down, table = "validated")

# Remove redundant target genes
target_results_up_unique <- unique(target_results_up@data$target_entrez)
target_results_down_unique <- unique(target_results_down@data$target_entrez)


for (mirna_symbol in dose24_miRNA_list_up) {
  
  targets <- get_multimir(mirna = mirna_symbol, table = 'mirtarbase')
  
  genes <- targets@data$target_entrez
  
  genes_unique <- unique((targets@data$target_entrez))
  
  print(mirna_symbol)
  
  print('Genes:')
  
  print(length(genes_unique))
  
  ekegg <- enrichKEGG(genes_unique)
  
  print('Pathways:')
  
  print(length(ekegg@result$p.adjust[ekegg@result$p.adjust < 1]))
  
  print('Pathways below Padj 0.05:')
  
  print(length(ekegg@result$p.adjust[ekegg@result$p.adjust < 0.05]))
#  result <- as.data.frame(ekegg)
  
#  print(result)
 # return(result)

}

tgtscn <- get_multimir(mirna = 'hsa-miR-26a-5p', table = 'targetscan')

mirdb <- get_multimir(mirna = 'hsa-miR-26a-5p', table = 'mirdb')
length(unique(tgtscn@data$target_symbol))

mirtarbase <- get_multimir(mirna = 'hsa-miR-26a-5p', table = 'mirtarbase')
length(unique(mirtarbase@data$target_symbol))

tarbase <- get_multimir(mirna = 'hsa-miR-26a-5p', table = 'tarbase')
length(unique(tarbase@data$target_symbol))

length(unique(mirdb@data$target_symbol))
View(mirdb@data)

table(targets@data$database)

targets@data %>%
  filter(database == 'tarbase') %>%
  head(20)

targets@data %>%
  filter(database == 'mirtarbase') %>%
  head(20)

targets@data %>%
  filter(database == 'mirecords') %>%
  head(20)

#==========================#
# KEGG Enrichment
#==========================#


# Isolate target genes present in airway background universe 

# Upregulated
genes_entrez_up <- intersect(
  unique(target_results_up_filtered$target_entrez),
  hbec_genes_entrez
)

# Downregulated
genes_entrez_down <- intersect(
  unique(target_results_down_filtered$target_entrez),
  hbec_genes_entrez
)

# Perform KEGG enrichment analysis
ekegg_up <- enrichKEGG(genes_entrez_up, universe = hbec_genes_entrez)
ekegg_down <- enrichKEGG(genes_entrez_down, universe = hbec_genes_entrez)

# Visualise KEGG enrichment analysis results
dotplot(ekegg_up)
dotplot(ekegg_down)

ekegg_up@result %>% 
  dplyr::select(p.adjust, Description, category, subcategory) %>%
  head()

ekegg_up@result %>% View()






