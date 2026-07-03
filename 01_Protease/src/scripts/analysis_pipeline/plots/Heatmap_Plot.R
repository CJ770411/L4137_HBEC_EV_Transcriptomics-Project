# Load necessary libraries
library(edgeR)
library(dplyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)


#==========================#
# Load DE object     
#==========================#


all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose24.rds"
)

dge <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/dge.rds"
)


#==========================#
# Heatmap Plot 
#==========================#

# Create subset of significant DE miRNAs 
top_mirnas <- all_results_dose24[all_results_dose24$FDR<0.05,]

# Calculate log2CPM values
logCPM <- cpm(dge, log = TRUE, prior.count = 1) # Adds 1 count to all observations to prevent inf values

# Set miRNA symbols as row names
rownames(logCPM) <- dge$genes$Symbol

# Refine to only subset of top miRNAs
logCPM_subset <- logCPM[top_mirnas$Symbol, ]

# Scale for z-score (relative expression levels)
logCPM_subset_scaled <- t(scale(t(logCPM_subset)))

# Extract groups (treatment conditions) from 'dge' object
group <- dge$samples$group

# Create annotation column containing groupings
annotation_col <- data.frame(Group = group)

# Set sample ID as row names
rownames(annotation_col) <- colnames(logCPM_subset_scaled)

# Create colour palette
heatmap_colours <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
  )(100)

# Create heatmap
pheatmap(
  logCPM_subset_scaled,
  annotation_col = annotation_col,
  scale = 'none',
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  show_rownames = TRUE,
  fontsize_col = 13,
  border_color = NA,
  color = heatmap_colours
)


