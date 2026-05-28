library(readr)
library(dplyr)

# Read in results from '04b_DE_analysis_edgeR'
dose24_sig_data <- read_csv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv")
count_matrix <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/deseq2_count_matrix.csv")

# Remove rows with 0 counts
count_matrix_trimmed <- count_matrix[rowSums(count_matrix[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors

# Set miRNA symbol column name for clarity
colnames(count_matrix_trimmed)[1] <- "Symbol"

# Sort significant miRNAs based on logFC
dose24_sig_mirna_ranked_lfc_all <- dose24_sig_data %>% dplyr::arrange(desc(abs(logFC)))

# Extract UP- and DOWN-regulated miRNAs
dose24_sig_mirna_ranked_lfc_up <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC > 0) # Upregulated
dose24_sig_mirna_ranked_lfc_down <- dose24_sig_mirna_ranked_lfc_all %>% filter(logFC < 0) # Downregulated


# Create miRNA lists
dose24_sig_mirnas_ranked_all <- dose24_sig_mirna_ranked_lfc_all$Symbol # All significant miRNAs
dose24_sig_mirnas_ranked_up <- dose24_sig_mirna_ranked_lfc_up$Symbol # Upregulated miRNAs
dose24_sig_mirnas_ranked_down <- dose24_sig_mirna_ranked_lfc_down$Symbol # Downregulated miRNAs
dose24_background_mirnas <- count_matrix_trimmed$Symbol # All miRNAs with at least one read in one sample

# Save miRNA lists as TXT file
write(dose24_sig_mirnas_ranked_all, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/dose24_sig_mirnas_ranked_all.txt")
write(dose24_sig_mirnas_ranked_up, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/dose24_sig_mirnas_ranked_up.txt")
write(dose24_sig_mirnas_ranked_down, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/dose24_sig_mirnas_ranked_down.txt")
write(dose24_background_mirnas, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/dose24_background_mirnas.txt")









