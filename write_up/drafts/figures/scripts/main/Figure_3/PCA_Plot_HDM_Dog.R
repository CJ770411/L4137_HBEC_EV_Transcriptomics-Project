# Load necessary libraries
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)


#==========================#
# Load DE objects      
#==========================#


dge <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/dge.rds"
)

group <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/dge_group.rds"
)


#==========================#
# PCA Plot      
#==========================#


# Log2 CPM transformation 
logCPM <- cpm(dge, log = TRUE, prior.count = 2)

# Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE) # Transpose matrix 

# Calculate variance explained
var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# Create data frame for plotting
pca_data <- data.frame(
  Sample = colnames(logCPM),
  Group = group,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)


# Define PCA plot point colours
pca_cols <- c(
  "#E69F00",
  "#56B4E9", 
  "#009E73", 
  "#D55E00"
)


# Create PCA plot labels 
pca_data <- pca_data %>%
  mutate(Sample_label = Sample,
         Sample_label = str_replace_all(Sample_label, "__1", " (1)"),
         Sample_label = str_replace_all(Sample_label, "__2", " (2)"),
         Sample_label = str_replace_all(Sample_label, "__3", " (3)"),
         Sample_label = str_replace_all(Sample_label, "_", " "),
         Sample_label = str_replace_all(Sample_label, "Unstimulated", "Untreated")
         )


# Create PCA plot
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=Group, label = Sample_label)) +
  geom_point(size=9) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) + 
  scale_color_manual(labels = c('Unstimulated' = 'Untreated',
                                'Der_p1_Active' = 'Der p1 Active',
                                'Der_p1_Inactive' = 'Der p1 Inactive',
                                'Can_F_1' = 'Can f 1'),
                     values = pca_cols,
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
print(pca_plot)

# Save PCA plot
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/panels/PCA_plot_hdm_dog.svg", plot = pca_plot, width=10, height=7, units='in')




