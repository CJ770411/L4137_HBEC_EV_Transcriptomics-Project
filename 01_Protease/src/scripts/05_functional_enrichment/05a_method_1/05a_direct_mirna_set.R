library(readr)

# Read in results from '04b_DE_analysis_edgeR'
dose24_sig_data <- read_csv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv")
count_matrix <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/deseq2_count_matrix.csv")


# Remove rows with 0 counts
reads <- count_matrix[rowSums(count_matrix[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors


# [start from here]




### Generate miRNA lists ###
dose24_sig_data_ordered_lfc <- dose24_sig_data %>% dplyr::arrange(desc(abs(logFC)))
dose24_all_data_ordered_lfc <- dose24_all_data %>% dplyr::arrange(desc(abs(logFC)))

## Upregulated

miRNA_ranked_list_dose24_sig <- dose24_sig_data_ordered_lfc$Symbol
miRNA_ranked_list_dose24_all <- dose24_all_data_ordered_lfc$Symbol


write(miRNA_ranked_list_dose24_sig, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/miRNA_ranked_list_dose24_sig.txt")
write(miRNA_ranked_list_dose24_all, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/05_functional_enrichment/05a_method_1/miRNA_ranked_list_dose24_all.txt")






