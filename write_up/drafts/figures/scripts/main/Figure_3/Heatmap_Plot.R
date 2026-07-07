# Load necessary libraries
library(edgeR)
library(dplyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)


#==========================#
# Load DE object     
#==========================#


all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/all_results_dose24.rds"
)

dge <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/dge.rds"
)

#==========================#
# Heatmap Plot Data Prep
#==========================#


# Create subset of significant DE miRNAs 
top_mirnas <- all_results_dose24[all_results_dose24$FDR<0.05 & abs(all_results_dose24$logFC) > 1,]

# Calculate log2CPM values
logCPM <- cpm(dge, log = TRUE, prior.count = 1) # Adds 1 count to all observations to prevent inf values

# Set miRNA symbols as row names
rownames(logCPM) <- dge$genes$Symbol

# Refine to only subset of top miRNAs
logCPM_subset <- logCPM[top_mirnas$Symbol, ]

# Extract groups (treatment conditions) from 'dge' object
group <- dge$samples$group

# Create annotation column containing groupings
annotation_col <- data.frame(Group = group)

# Set sample ID as annotation column row names 
rownames(annotation_col) <- colnames(logCPM_subset)

# Create colour palette
heatmap_colours <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)


#==========================#
# Heatmap Plot (Dose24 vs Control)
#==========================#

logCPM_subset_d24 <- logCPM_subset[, c(
  "Control__1", "Control__2", "Control__3",
  "Dose24__1", "Dose24__2", "Dose24__3"
)]

logCPM_subset_d24_scaled <- t(scale(t(logCPM_subset_d24)))


# Create heatmap of 'Untreated' vs 'Dose 24'
untreated_vs_dose24_heatmap <- pheatmap(
  logCPM_subset_d24_scaled,
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

# Visualise heatmap
print(untreated_vs_dose24_heatmap)

# Save heatmap 
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/panels/untreated_vs_dose24_heatmap.svg", plot = untreated_vs_dose24_heatmap, width=10, height=7, units='in')


#==========================#
# Heatmap Plot (All conditions)
#==========================#

# Scale for z-score (relative expression levels)
logCPM_subset_scaled <- t(scale(t(logCPM_subset)))

# Create heatmap showing all conditions
all_contrasts_heatmap <- pheatmap(
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

# Visualise heatmap
print(all_contrasts_heatmap)

# Save heatmap
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/panels/all_contrasts_heatmap.svg", plot = all_contrasts_heatmap, width=10, height=7, units='in')


