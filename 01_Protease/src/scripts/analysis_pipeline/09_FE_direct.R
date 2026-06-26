# Load necessary libraries
library(readr)
library(dplyr)

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold

#==========================#
# Create Background Universe 
#==========================#

# Read in count matrix containing all contrasts
count_matrix <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/06_create_count_matrix/edgeR_count_matrix.txt")

# Isolate reads from Dose 24 vs Control
dose24_data  <- count_matrix %>% dplyr::select(Symbol,
                                               Control__1, Control__2, Control__3, 
                                               Dose24__1, Dose24__2, Dose24__3)


# Remove rows with 0 counts
dose24_data_filtered <- dose24_data[rowSums(dose24_data[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors

# Create list of miRNA IDs
dose24_background_mirnas <- dose24_data_filtered$Symbol

# Save background miRNA list as TXT file
write(dose24_background_mirnas, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_background_mirnas.txt")


#==========================#
# Create miRNA lists
#==========================#

# Read in results from '04b_DE_analysis_edgeR'
dose24_sig_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Sort significant miRNAs based on logFC
dose24_sig_mirna_ranked_lfc_all <- dose24_sig_data %>% dplyr::arrange(desc(abs(logFC)))

# Extract UP- and DOWN-regulated miRNAs
dose24_sig_mirna_ranked_lfc_up <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC > 0) # Upregulated
dose24_sig_mirna_ranked_lfc_down <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC < 0) # Downregulated

# Exclude miRNAs below logFC threshold
# Note: the 'all' list is used for (G)SEA where a log threshold isn't required as it uses the order
dose24_sig_mirna_ranked_lfc_up_filtered <- dose24_sig_mirna_ranked_lfc_up[dose24_sig_mirna_ranked_lfc_up$logFC>=logFC_threshold,]
dose24_sig_mirna_ranked_lfc_down_filtered <- dose24_sig_mirna_ranked_lfc_down[dose24_sig_mirna_ranked_lfc_down$logFC<=-logFC_threshold,]


# Create miRNA lists
dose24_sig_mirnas_ranked_all <- dose24_sig_mirna_ranked_lfc_all_filtered$Symbol # All significant miRNAs
dose24_sig_mirnas_ranked_up <- dose24_sig_mirna_ranked_lfc_up_filtered$Symbol # Upregulated miRNAs
dose24_sig_mirnas_ranked_down <- dose24_sig_mirna_ranked_lfc_down_filtered$Symbol # Downregulated miRNAs

# Save miRNA lists as TXT file
write(dose24_sig_mirnas_ranked_all, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_sig_mirnas_ranked_all.txt")
write(dose24_sig_mirnas_ranked_up, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_sig_mirnas_ranked_up.txt")
write(dose24_sig_mirnas_ranked_down, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_sig_mirnas_ranked_down.txt")









