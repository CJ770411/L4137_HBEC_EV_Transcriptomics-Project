# Load necessary libraries
library(readr)
library(dplyr)

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold

#==========================#
# Create miRNA lists
#==========================#

# Read in results from '04b_DE_analysis_edgeR'
dose24_sig_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv")
dose24_background_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose24.csv")

# Sort significant miRNAs based on logFC
dose24_sig_mirna_ranked_lfc_all <- dose24_sig_data %>% dplyr::arrange(desc(abs(logFC)))

# Extract UP- and DOWN-regulated miRNAs
dose24_sig_mirna_ranked_lfc_up <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC > 0) # Upregulated
dose24_sig_mirna_ranked_lfc_down <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC < 0) # Downregulated

# Exclude miRNAs below logFC threshold
dose24_sig_mirna_ranked_lfc_all_filtered <- dose24_sig_mirna_ranked_lfc_all[abs(dose24_sig_mirna_ranked_lfc_all$logFC)>=logFC_threshold,]
dose24_sig_mirna_ranked_lfc_up_filtered <- dose24_sig_mirna_ranked_lfc_up[dose24_sig_mirna_ranked_lfc_up$logFC>=logFC_threshold,]
dose24_sig_mirna_ranked_lfc_down_filtered <- dose24_sig_mirna_ranked_lfc_down[dose24_sig_mirna_ranked_lfc_down$logFC<=-logFC_threshold,]


# Create miRNA lists
dose24_sig_mirnas_ranked_all <- dose24_sig_mirna_ranked_lfc_all_filtered$Symbol # All significant miRNAs
dose24_sig_mirnas_ranked_up <- dose24_sig_mirna_ranked_lfc_up_filtered$Symbol # Upregulated miRNAs
dose24_sig_mirnas_ranked_down <- dose24_sig_mirna_ranked_lfc_down_filtered$Symbol # Downregulated miRNAs
dose24_background_mirnas <- dose24_background_data$Symbol # All miRNAs with at least one read in one sample

# Save miRNA lists as TXT file
write(dose24_sig_mirnas_ranked_all, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/09_FE_direct/dose24_sig_mirnas_ranked_all.txt")
write(dose24_sig_mirnas_ranked_up, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/09_FE_direct/dose24_sig_mirnas_ranked_up.txt")
write(dose24_sig_mirnas_ranked_down, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/09_FE_direct/dose24_sig_mirnas_ranked_down.txt")
write(dose24_background_mirnas, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/09_FE_direct/dose24_background_mirnas.txt")









