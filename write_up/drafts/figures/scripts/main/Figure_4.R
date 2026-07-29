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
# Panel C - Differential Expression Barchart
#==========================#


# Generate Bar Chart Data (Protease)   

# Load primary HBEC DE data
dose0_005_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose0005.csv")
dose0_08_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose008.csv")
dose2_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose2.csv")
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Add condition labels
dose0_005_data <- dose0_005_data %>% 
  count() %>%
  mutate(Condition = 'Dose 0.005')

dose0_08_data <- dose0_08_data %>% 
  count() %>%
  mutate(Condition = 'Dose 0.08')

dose2_data <- dose2_data %>% 
  count() %>%
  mutate(Condition = 'Dose 2')

dose24_data <- dose24_data %>% 
  count() %>%
  mutate(Condition = 'Dose 24')



# Generate Bar Chart Data (HDM_Dog)   

# Load DE data
can_f_1_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Can_F_1.csv")
der_p1_active_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Active.csv")
der_p1_inactive_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Inactive.csv")

# Add condition labels
can_f_1_data <- can_f_1_data %>% 
  count() %>%
  mutate(Condition = 'Can f 1')

der_p1_active_data <- der_p1_active_data %>% 
  count() %>%
  mutate(Condition = 'Der P1 Active')

der_p1_inactive_data <- der_p1_inactive_data %>% 
  count() %>%
  mutate(Condition = 'Der P1 Inactive')



# Create Bar Chart 

# Combine condition data
boxplot_data <- rbind(dose0_005_data,
                      dose0_08_data,
                      dose2_data,
                      dose24_data,
                      can_f_1_data,
                      der_p1_active_data,
                      der_p1_inactive_data)

# Preserve condition order
boxplot_data$Condition <- factor(
  boxplot_data$Condition,
  levels = c(
    "Dose 0.005",
    "Dose 0.08",
    "Dose 2",
    "Dose 24",
    "Can f 1",
    "Der P1 Active",
    "Der P1 Inactive"
  )
)

# Set column name as DE miRNA count
colnames(boxplot_data)[colnames(boxplot_data) == 'n'] <- 'Count'

# Plot bar chart
diffexpr_barchart <- ggplot(boxplot_data, aes(x = Condition, y = Count, fill = Condition)) +
  geom_col() +
  labs(x = 'Condition',
       y = 'DE miRNA Count') +
  theme_classic() + 
  labs_pubr() +
  scale_fill_manual(values = c(
    'Dose 0.005'      = '#D1EEE8',
    'Dose 0.08'       = '#7BC8B6',
    'Dose 2'          = '#2A9D8F',
    'Dose 24'         = '#146C63',
    'Can f 1'         = '#E8D5F2',
    'Der P1 Active'   = '#C7A0D8',
    'Der P1 Inactive' = '#9C6BB3',
    'Poly I:C'        = '#6A3D9A'
  )) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_y_continuous(
    breaks = seq(0, 30, by = 5))

# Visualise bar chart
print(diffexpr_barchart)



#==========================#
# Combine Plot to Create Figure
#==========================#

# Combine plots
figure_4 <- ((primaryHBEC_pca_plot | calu3_pca_plot) /  diffexpr_barchart) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))

# Visualise figure
figure_4

# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_4/figure_4.pdf", plot = figure_4, height=7, width=10, dpi=600, units="in", scale=2)


