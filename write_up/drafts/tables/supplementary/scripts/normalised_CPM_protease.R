# Load necessary libraries
library(edgeR)

# This produces normalised CPM values for all miRNAs included in the differential expression analysis.
# This means all the miRNAs retained after filtering.

# Load edgeR DE object 
dge_protease <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/dge.rds"
)

# Calculate CPM-normalised counts
CPM_protease <- cpm(dge_protease, log = FALSE, prior.count = 1) # Adds 1 count to all observations to prevent inf values

# Set miRNA symbols as row names
rownames(CPM_protease) <- dge_protease$genes$Symbol

# Save normalised counts as CSV
write.csv(CPM_protease, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/tables/supplementary/normalised_CPM_Protease.csv")
