# Load necessary libraries
library(readr)
library(dplyr)

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold for ORA
FDR_threshold <- 0.05 # FDR significance threshold for ORA

#==========================#
# Read in edgeR data
#==========================#

# Read in all DE miRNAs from edgeR
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/07_DE_analysis/all_results_dose24.csv")


#==========================#
# Create Background Universe 
#==========================#

# Create list of miRNA IDs
dose24_background_mirnas <- dose24_data$Symbol

# Save background miRNA list as TXT file
write(dose24_background_mirnas, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_background.txt")


#==========================#
# Create miRNA lists
#==========================#

### ORA  ###

# Extract significant UP- and DOWN-regulated miRNAs
dose24_ora_up <- dose24_data %>% filter(logFC >= logFC_threshold, FDR < FDR_threshold) %>% pull(Symbol) # Upregulated
dose24_ora_down <- dose24_data %>% filter(logFC <= -logFC_threshold, FDR < FDR_threshold) %>% pull(Symbol) # Downregulated
dose24_ora_all <- dose24_data %>% filter(abs(logFC) >= logFC_threshold, FDR < FDR_threshold) %>% pull(Symbol) # All

# Save miRNA lists as TXT file
write(dose24_ora_up, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_ora_up.txt")
write(dose24_ora_down, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_ora_down.txt")
write(dose24_ora_all, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_ora_all.txt")



### (G)SEA  ###
# Note: the 'all' list is used for (G)SEA where a log threshold isn't required as it uses the order of logFC expression

# Sort significant miRNAs based on logFC
dose24_gsea <- dose24_data %>% dplyr::arrange(desc(logFC)) %>% pull(Symbol)

# Save miRNA lists as TXT file
write(dose24_gsea, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/09_FE_direct/dose24_gsea.txt")




