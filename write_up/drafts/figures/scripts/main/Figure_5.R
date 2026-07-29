# Load necessary libraries
library(readr)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplotify)
library(patchwork)


#==========================#
# Panel A - MA Plot
#==========================#

# User Configuration 
logFC_threshold <- 1 # miRNA expression threshold
FDR_threshold <- 0.05 # Significance threshold


# Load required data  
# All miRNAs tested for differential expression in primary HBEC dose 24ug condition
all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/all_results_dose24.rds"
)


# Determine significance to allow colour-by-significance in the plot
all_results_dose24$significance <- "Not significant"
all_results_dose24$significance[all_results_dose24$FDR < FDR_threshold & all_results_dose24$logFC > logFC_threshold] <- "Upregulated"
all_results_dose24$significance[all_results_dose24$FDR < FDR_threshold & all_results_dose24$logFC < -logFC_threshold] <- "Downregulated"

# Preserve the order of appearance in legend of MA plot
all_results_dose24$significance <- factor(all_results_dose24$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

# Set axis limits
x_min <- floor(min(all_results_dose24$logCPM))
x_max <- ceiling(max(all_results_dose24$logCPM))

# Create MA plot
ma_plot <- ggplot(all_results_dose24, aes(x = logCPM, y = logFC, color = significance, size = significance, alpha = significance)) +
  geom_point() +
  scale_color_manual(
    values = c("Not significant" = "grey80",
               "Downregulated" = "#0072B2",
               "Upregulated" = "#E64B35")
  ) +
  scale_size_manual(    
    values = c("Not significant" = 3,
               "Downregulated" = 4,
               "Upregulated" = 4)) +
  scale_alpha_manual(    
    values = c("Not significant" = 0.4,
               "Downregulated" = 0.9,
               "Upregulated" = 0.9)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  labs(x = "Average log2 CPM",
       y = "Log2 Fold Change",
       color = '',
       size = '',
       alpha = '') +
  theme_pubclean() +
  labs_pubr() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(ylim = c(10, -10)) +
  scale_x_continuous(
    breaks = seq(x_min, x_max, by = 1)) +
  scale_y_continuous(
    breaks = seq(-7.5, 7.5, by = 2.5)
  )

# Visualise MA plot
print(ma_plot)



#==========================#
# Panel B - Volcano Plot
#==========================#

# User Configuration 
logFC_threshold <- 1 # miRNA expression threshold
FDR_threshold <- 0.05 # Significance threshold


# Load required data  
# All miRNAs tested for differential expression in primary HBEC dose 24ug condition  
all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/all_results_dose24.rds"
)

# Determine significance to allow colour-by-significance in the plot
all_results_dose24$significance <- "Not significant"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC > logFC_threshold] <- "Upregulated"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC < -logFC_threshold] <- "Downregulated"

# Preserve the order of appearance in legend of MA plot
all_results_dose24$significance <- factor(all_results_dose24$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

# Set axis limits
x_max <- ceiling(max(abs(all_results_dose24$logFC)))
x_min <- -x_max

# Create volcano plot
volcano_plot <- ggplot(all_results_dose24, aes(x = logFC, y = -log10(FDR), color = significance, alpha = significance, size = significance)) +
  geom_point() +
  theme_classic() +
  labs_pubr() +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = c(1, -1), col = 'black', linetype = 'dashed', alpha = 0.6) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white')) +
  labs(x = "Log2 Fold Change", 
       y = "-log10 FDR",
       color = '',
       alpha = '',
       size = '') +
  scale_color_manual(values = c('Not significant' = 'grey80', 'Upregulated' = '#E64B35', 'Downregulated' = '#0072B2')) +
  scale_size_manual(    
    values = c("Not significant" = 3,
               "Downregulated" = 4,
               "Upregulated" = 4)) +
  scale_alpha_manual(    
    values = c("Not significant" = 0.4,
               "Downregulated" = 0.9,
               "Upregulated" = 0.9)) +
  coord_cartesian(xlim = c(x_max, x_min)) +
  scale_x_continuous(
    breaks = seq(-10, 10, by = 2.5)
  ) 

# Visualise volcano plot
print(volcano_plot)




#==========================#
# Panel C - Heatmap
#==========================#

# Load required data   
# DGE object
primaryHBEC_dge <- readRDS("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/dge.rds")

# Significant DE miRNAs from primary HEBC dose 24ug exposure
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Calculate log2CPM values
heatmap_logCPM <- cpm(primaryHBEC_dge, log = TRUE, prior.count = 1) # Adds 1 count to all observations to prevent inf values

# Set miRNA symbols as row names
rownames(heatmap_logCPM) <- primaryHBEC_dge$genes$Symbol

# Refine to only subset of top miRNAs
heatmap_logCPM_subset <- heatmap_logCPM[dose24_data$Symbol, ]

# Extract groups (treatment conditions) from dge object
heatmap_group <- primaryHBEC_dge$samples$group

# Create annotation column
heatmap_annotation_col <- data.frame(Condition = c('Untreated', 'Untreated', 'Untreated',
                                                   'Dose 24', 'Dose 24', 'Dose 24'))

# Set sample ID as annotation column row names 
rownames(heatmap_annotation_col) <- c('Untreated (1)', 'Untreated (2)', 'Untreated (3)',
                                      'Dose 24 (1)', 'Dose 24 (2)', 'Dose 24 (3)')

# Create subset containing only samples in the dose 24 experiment
heatmap_logCPM_subset_d24 <- heatmap_logCPM_subset[, c(
  "Control__1", "Control__2", "Control__3",
  "Dose24__1", "Dose24__2", "Dose24__3"
)]

colnames(heatmap_logCPM_subset_d24) <- c('Untreated (1)', 'Untreated (2)', 'Untreated (3)',
                                         'Dose 24 (1)', 'Dose 24 (2)', 'Dose 24 (3)')

# Transpose and standardise expression across samples
logCPM_subset_d24_scaled <- t(scale(t(heatmap_logCPM_subset_d24)))

# Create colour palette
heatmap_colours <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

# Create heatmap of 'Untreated' vs 'Dose 24'
untreated_vs_dose24_heatmap <- pheatmap(
  logCPM_subset_d24_scaled,
  annotation_col = heatmap_annotation_col,
  scale = 'none',
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  show_rownames = TRUE,
  fontsize_col = 13,
  border_color = NA,
  color = heatmap_colours
)

# Convert heatmap to ggplot object for compatibility with patchwork functions
untreated_vs_dose24_heatmap_ggplot <- as.ggplot(untreated_vs_dose24_heatmap)

# Visualise heatmap
print(untreated_vs_dose24_heatmap_ggplot)



#==========================#
# Combine Plot to Create Figure
#==========================#

# Combine plots
figure_5 <- ((ma_plot / volcano_plot) | untreated_vs_dose24_heatmap_ggplot) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))

# Visualise figure
print(figure_5)

# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_5/figure_5.pdf", plot = figure_5, height=8, width=12, dpi=600, units="in", scale=2)


