# Load necessary libraries
library(readr)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)


### OBSOLETE 18.06.2026 due to removal of airway-specific filtering and combined pathway analysis. New method to be adopted is individual (per-miRNA) pathway analysis.

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1.5 # miRNA expression threshold
mirna_count <- 3 # Minimum number of miRNAs required to target a target gene
database_count <- 1 # Minimum number of databases a target gene must appear in


#==========================#
# Create Background Universe 
#==========================#


# Gene expression data from https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data#cell_type_data 
# download file: rna_single_cell_type_cell_types.tsv.zip
hbec_gene_file <- read.delim("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/rna_single_cell_type.tsv")

# Create list of HBEC-related cells
airway_cells <- c(
  "respiratory basal cells",
  "respiratory ciliated cells",
  "respiratory secretory cells",
  "respiratory ionocytes",
  "respiratory deuterosomal cells"
)

# Get unique genes from HBEC-related cell types and have >1 count per million reads
hbec_genes <- hbec_gene_file %>%
  filter(`Cell.type` %in% airway_cells,
         nCPM > 1) %>%
  pull(`Gene.name`) %>%
  unique()

# Convert gene symbol to Entrez format
hbec_genes_entrez <- bitr(
  hbec_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Convert Entrez genes to character vector for KEGG analysis
hbec_genes_entrez <- hbec_genes_entrez$ENTREZID



#==========================#
# Identify Target Genes
#==========================#

# Read in signficant results from dose 24 vs control
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Extract UP- and DOWN-regulated miRNAs
dose24_data_up <- dose24_data %>% filter(logFC > logFC_threshold) # Upregulated
dose24_data_down <- dose24_data %>% filter(logFC < -logFC_threshold) # Downregulated

# Create list of miRNA IDs
dose24_miRNA_list_up <- dose24_data_up$Symbol
dose24_miRNA_list_down <- dose24_data_down$Symbol

# Identify gene targets for each miRNA
target_results_up <- get_multimir(mirna = dose24_miRNA_list_up, table = "validated")
target_results_down <- get_multimir(mirna = dose24_miRNA_list_down, table = "validated")

# Filter target genes based on:
#     1. No. miRNAs targeting the target gene
#     2. No. databases the target gene appears in

# Upregulated
target_results_up_filtered <- target_results_up@data %>%
    group_by(target_entrez) %>%
    summarise(
      n_miRNAs = n_distinct(mature_mirna_id),
      n_db = n_distinct(database)
    ) %>%
    filter(n_miRNAs >= mirna_count,
           n_db >= database_count)

# Downregulated
target_results_down_filtered <- target_results_down@data %>%
  group_by(target_entrez) %>%
  summarise(
    n_miRNAs = n_distinct(mature_mirna_id),
    n_db = n_distinct(database)
  ) %>%
  filter(n_miRNAs >= mirna_count,
         n_db >= database_count)


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






