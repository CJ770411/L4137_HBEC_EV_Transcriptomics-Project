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

#==========================#
# Create Background Universe 
#==========================#

# Read in count matrix containing all contrasts
count_matrix <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/06_create_count_matrix/edgeR_count_matrix.txt")

# Isolate reads from Dose 24 vs Control
dose24_data  <- count_matrix %>% dplyr::select(Symbol,
                                        Control__1, Control__2, Control__3, 
                                        Dose24__1, Dose24__2, Dose24__3)


# Remove rows with 0 counts
dose24_data_filtered <- dose24_data[rowSums(dose24_data[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors

# Create list of miRNA IDs
dose24_mirnas_all <- dose24_data_filtered$Symbol

# Identify gene targets for all DE miRNAs
target_results_background <- get_multimir(mirna = dose24_mirnas_all, table = "mirtarbase")

# Remove redundant target genes
target_results_background_unique <- unique(target_results_background@data$target_entrez) # Entrez format genes


#==========================#
# Identify Target Genes
#==========================#

# Read in signficant results from dose 24 vs control
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv")

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

# Perform pathway enrichment
epath_up <- enrichPathway(target_results_up_unique, universe = target_results_background_unique)
epath_down <- enrichPathway(target_results_down_unique, universe = target_results_background_unique)

# Visualise pathway enrichment
dotplot(epath_up)
dotplot(epath_down)

# View pathway enrichment as dataframes for exploration
#epath_up@result %>% View()
#epath_down@result %>% View()







