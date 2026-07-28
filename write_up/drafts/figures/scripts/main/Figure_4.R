# Load necessary libraries
library(readr)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(patchwork)



#==========================#
# Panel A - PCA (Primary HBECs)
#==========================#

# Load required DE objects      
primaryHBEC_dge <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/dge.rds"
)

primaryHBEC_group <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/dge_group.rds"
)


# Perform log2 CPM transformation 
primaryHBEC_logCPM <- cpm(primaryHBEC_dge, log = TRUE, prior.count = 2)

# Run PCA
primaryHBEC_pca <- prcomp(t(primaryHBEC_logCPM), scale. = TRUE) # Transpose matrix 

# Calculate variance explained
primaryHBEC_var_explained <- (primaryHBEC_pca$sdev^2) / sum(primaryHBEC_pca$sdev^2) * 100

# Create data frame for plotting
primaryHBEC_pca_data <- data.frame(
  Sample = colnames(primaryHBEC_logCPM),
  Group = primaryHBEC_group,
  PC1 = primaryHBEC_pca$x[,1],
  PC2 = primaryHBEC_pca$x[,2]
)

# Create PCA plot labels 
primaryHBEC_pca_data <- primaryHBEC_pca_data %>%
  mutate(Sample_label = Sample,
         Sample_label = str_replace(Sample_label, "Dose", ""),
         Sample_label = str_replace_all(Sample_label, "0_", "0."),
         Sample_label = str_replace_all(Sample_label, "_", " "),
         Sample_label = str_replace_all(Sample_label, " 1", "(1)"),
         Sample_label = str_replace_all(Sample_label, " 2", "(2)"),
         Sample_label = str_replace_all(Sample_label, " 3", "(3)"),
         Sample_label = str_replace_all(Sample_label, "Control", "Untreated")
  )


# Define PCA plot point colours
primaryHBEC_pca_cols <- c(
  "#E69F00",
  "#56B4E9", 
  "#009E73", 
  "#D55E00", 
  "#CC79A7"  
)


# Create PCA plot
primaryHBEC_pca_plot <- ggplot(primaryHBEC_pca_data, aes(PC1, PC2, color=Group, label = Sample_label)) +
  geom_point(size=9) +
  xlab(paste0("PC1 (", round(primaryHBEC_var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(primaryHBEC_var_explained[2], 1), "%)")) + 
  scale_color_manual(labels = c('Control' = 'Untreated',
                                'Dose_0_005' = 'Dose 0.005',
                                'Dose_0_08' = 'Dose 0.08',
                                'Dose_2' = 'Dose 2',
                                'Dose_24' = 'Dose 24'),
                     values = primaryHBEC_pca_cols,
                     name = 'Condition') +
  geom_text_repel(box.padding = 0.9) +
  theme_bw() +
  labs_pubr() +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

# Visualise PCA plot
print(primaryHBEC_pca_plot)


#==========================#
# Panel B - PCA (Calu-3)
#==========================#

# Load DE objects      
calu3_dge <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/dge.rds"
)

calu3_group <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/dge_group.rds"
)

# Perform log2 CPM transformation 
calu3_logCPM <- cpm(calu3_dge, log = TRUE, prior.count = 2)

# Run PCA
calu3_pca <- prcomp(t(calu3_logCPM), scale. = TRUE) # Transpose matrix 

# Calculate variance explained
calu3_var_explained <- (calu3_pca$sdev^2) / sum(calu3_pca$sdev^2) * 100

# Create data frame for plotting
calu3_pca_data <- data.frame(
  Sample = colnames(calu3_logCPM),
  Group = calu3_group,
  PC1 = calu3_pca$x[,1],
  PC2 = calu3_pca$x[,2]
)


# Define PCA plot point colours
calu3_pca_cols <- c(
  "#E69F00",
  "#56B4E9", 
  "#009E73", 
  "#D55E00"
)


# Create PCA plot labels 
calu3_pca_data <- calu3_pca_data %>%
  mutate(Sample_label = Sample,
         Sample_label = str_replace_all(Sample_label, "__1", " (1)"),
         Sample_label = str_replace_all(Sample_label, "__2", " (2)"),
         Sample_label = str_replace_all(Sample_label, "__3", " (3)"),
         Sample_label = str_replace_all(Sample_label, "_", " "),
         Sample_label = str_replace_all(Sample_label, "Unstimulated", "Untreated")
  )


# Create PCA plot
calu3_pca_plot <- ggplot(calu3_pca_data, aes(PC1, PC2, color=Group, label = Sample_label)) +
  geom_point(size=9) +
  xlab(paste0("PC1 (", round(calu3_var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(calu3_var_explained[2], 1), "%)")) + 
  scale_color_manual(labels = c('Unstimulated' = 'Untreated',
                                'Der_p1_Active' = 'Der p1 Active',
                                'Der_p1_Inactive' = 'Der p1 Inactive',
                                'Can_F_1' = 'Can f 1'),
                     values = calu3_pca_cols,
                     name = 'Condition') +
  geom_text_repel(box.padding = 0.9) +
  theme_bw() +
  labs_pubr() +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

# Visualise PCA plot
print(calu3_pca_plot)



#==========================#
# Combine Plot to Create Figure
#==========================#

# Combine plots
figure_4 <- (primaryHBEC_pca_plot | calu3_pca_plot) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))

# Visualise figure
figure_4

# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_4/figure_4.pdf", plot = figure_4, height=4, width=9, dpi=600, units="in", scale=2)


