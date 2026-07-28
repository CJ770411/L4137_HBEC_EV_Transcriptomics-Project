# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(scales)
library(patchwork)
library(tidyr)

#==========================#
# Read-in Data
#==========================#

# Load MultiQC data file for miRTrace RNA categories plot
mirtrace_data_primaryHBECs <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_rna_categories_plot.txt")
mirtrace_data_Calu3 <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_rna_categories_plot.txt")



#==========================#
# Define Experimental Condition Groups (Primary HBECs)
#==========================#

# Untreated
primaryHBEC_untreated <- c(
  "NGS-110-028_S28_R1_001",
  "NGS-110-016_S16_R1_001",
  "NGS-110-007_S7_R1_001"
)

# Protease Dose 0.005
primaryHBEC_dose_0_005 <- c(
  "NGS-110-029_S29_R1_001",
  "NGS-110-017_S17_R1_001",
  "NGS-110-011_S11_R1_001"
)

# Protease Dose 0.08
primaryHBEC_dose_0_08 <- c(
  "NGS-110-012_S12_R1_001",
  "NGS-110-025_S25_R1_001",
  "NGS-110-006_S6_R1_001"
)

# Protease Dose 2
primaryHBEC_dose_2 <- c(
  "NGS-110-004_S4_R1_001",
  "NGS-110-026_S26_R1_001",
  "NGS-110-003_S3_R1_001"
)

# Protease Dose 24
primaryHBEC_dose_24 <- c(
  "NGS-110-030_S30_R1_001",
  "NGS-110-019_S19_R1_001",
  "NGS-110-008_S8_R1_001"
)


#==========================#
# Define Experimental Condition Groups (Calu-3)
#==========================#

# Untreated
calu3_untreated <- c(
  "NGS-110-021_S21_R1_001",
  "NGS-110-018_S18_R1_001",
  "NGS-110-002_S2_R1_001"
)

# Der p 1 (HDM) - active
calu3_hdm_active <- c(
  "NGS-110-027_S27_R1_001",
  "NGS-110-005_S5_R1_001",
  "NGS-110-023_S23_R1_001"
)

# Der p 1 (HDM) - inactive
calu3_hdm_inactive <- c(
  "NGS-110-020_S20_R1_001",
  "NGS-110-009_S9_R1_001",
  "NGS-110-001_S1_R1_001"
)

# Can f 1 (dog)
calu3_dog <- c(
  "NGS-110-024_S24_R1_001",
  "NGS-110-015_S15_R1_001",
  "NGS-110-013_S13_R1_001"
)


#==========================#
# Panel A - Overall miRNA Abundance Boxplot
#==========================#

# Extract miRNA-assigned reads for all primary HBEC samples
mirna_reads_primaryHBECs <- mirtrace_data_primaryHBECs %>%
  mutate(Condition = 'Primary HBEC') %>%
  select(miRNA, Sample, Condition) 

# Extract miRNA-assigned reads for all primary HBEC samples
mirna_reads_Calu3 <- mirtrace_data_Calu3 %>%
  mutate(Condition = 'Calu-3') %>%
  select(miRNA, Sample, Condition) 

# Merge primary HBEC and Calu-3 data into one dataframe 
overall_mirna_boxplot_data <- rbind(mirna_reads_primaryHBECs,
                                    mirna_reads_Calu3)

# Create boxplot
overall_mirna_boxplot <- ggplot(overall_mirna_boxplot_data, aes(x = Condition, y = miRNA, fill = Condition)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  labs(x = 'Experiment',
       y = 'miRNA-assigned Read Count') +
  theme_pubr() +
  labs_pubr() +
  theme(axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.key.size = unit(0.7, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 150000)) +
  scale_y_continuous(
    labels = unit_format(scale = 1e-3, suffix = "K")
  ) 

# Visualise boxplot
print(overall_mirna_boxplot)




#==========================#
# Panel B - Per-condition miRNA Abundance Boxplot (Primary HBECs)
#==========================#

# Extract miRNA-assigned reads for each primary HBEC condition

# Untreated
mirna_reads_primaryHBECs_untreated <- mirtrace_data_primaryHBECs %>%
  filter(Sample %in% primaryHBEC_untreated) %>%
  mutate(Condition = 'Untreated') %>%
  select(miRNA, Sample, Condition)

# Dose 0.005
mirna_reads_primaryHBECs_0_005 <- mirtrace_data_primaryHBECs %>%
  filter(Sample %in% primaryHBEC_dose_0_005) %>%
  mutate(Condition = 'Dose 0.005') %>%
  select(miRNA, Sample, Condition)

# Dose 0.08
mirna_reads_primaryHBECs_0_08 <- mirtrace_data_primaryHBECs %>%
  filter(Sample %in% primaryHBEC_dose_0_08) %>%
  mutate(Condition = 'Dose 0.08') %>%
  select(miRNA, Sample, Condition)

# Dose 2
mirna_reads_primaryHBECs_2 <- mirtrace_data_primaryHBECs %>%
  filter(Sample %in% primaryHBEC_dose_2) %>%
  mutate(Condition = 'Dose 2') %>%
  select(miRNA, Sample, Condition)

# Dose 24
mirna_reads_primaryHBECs_24 <- mirtrace_data_primaryHBECs %>%
  filter(Sample %in% primaryHBEC_dose_24) %>%
  mutate(Condition = 'Dose 24') %>%
  select(miRNA, Sample, Condition)

# Merge primary HBEC conditions data into one dataframe 
primaryHBEC_mirna_boxplot_data <- rbind(mirna_reads_primaryHBECs_untreated,
                                        mirna_reads_primaryHBECs_0_005,
                                        mirna_reads_primaryHBECs_0_08,
                                        mirna_reads_primaryHBECs_2,
                                        mirna_reads_primaryHBECs_24)


# Create boxplot
primaryHBEC_mirna_boxplot <- ggplot(primaryHBEC_mirna_boxplot_data, aes(x = Condition, y = miRNA, fill = Condition)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  labs(x = 'Condition',
       y = 'miRNA-assigned Read Count') +
  theme_pubr() +
  labs_pubr() +
  theme(axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.key.size = unit(0.7, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 150000)) +
  scale_y_continuous(
    labels = unit_format(scale = 1e-3, suffix = "K")
  ) 

# Visualise boxplot
print(primaryHBEC_mirna_boxplot)



#==========================#
# Panel C - Per-condition miRNA Abundance Boxplot (Calu-3)
#==========================#

# Extract miRNA-assigned reads for each Calu-3 condition

# Untreated
mirna_reads_calu3_untreated <- mirtrace_data_Calu3 %>%
  filter(Sample %in% calu3_untreated) %>%
  mutate(Condition = 'Untreated') %>%
  select(miRNA, Sample, Condition)

# Der p 1 (HDM) - active
mirna_reads_calu3_hdm_active <- mirtrace_data_Calu3 %>%
  filter(Sample %in% calu3_hdm_active) %>%
  mutate(Condition = 'HDM Active') %>%
  select(miRNA, Sample, Condition)

# Der p 1 (HDM) - inactive
mirna_reads_calu3_hdm_inactive <- mirtrace_data_Calu3 %>%
  filter(Sample %in% calu3_hdm_inactive) %>%
  mutate(Condition = 'HDM Inactive') %>%
  select(miRNA, Sample, Condition)

# Can f 1 (dog)
mirna_reads_calu3_dog <- mirtrace_data_Calu3 %>%
  filter(Sample %in% calu3_dog) %>%
  mutate(Condition = 'Dog') %>%
  select(miRNA, Sample, Condition)

# Merge Calu-3 conditions data into one dataframe 
calu3_mirna_boxplot_data <- rbind(mirna_reads_calu3_untreated,
                              mirna_reads_calu3_hdm_active,
                              mirna_reads_calu3_hdm_inactive,
                              mirna_reads_calu3_dog)

# Create boxplot
calu3_mirna_boxplot <- ggplot(calu3_mirna_boxplot_data, aes(x = Condition, y = miRNA, fill = Condition)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  labs(x = 'Condition',
       y = 'miRNA-assigned Read Count') +
  theme_pubr() +
  labs_pubr() +
  theme(axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.key.size = unit(0.7, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 5000)) + 
  scale_y_continuous(
    labels = unit_format(scale = 1e-3, suffix = "K")
  )

# Visualise boxplot
print(calu3_mirna_boxplot)


#==========================#
# Panel D - RNA Composition Plot (Primary HBECs)
#==========================#

# Calculate library size to enable proportion calculations 
mirtrace_data_primaryHBECs <- mirtrace_data_primaryHBECs %>%
  mutate(Total = rowSums(mirtrace_data_primaryHBECs[-1]))

# Calculate proportion of each RNA type as percentages
mirtrace_data_primaryHBECs_percent <- mirtrace_data_primaryHBECs %>%
  mutate(
    miRNA = miRNA / Total * 100,
    rRNA = rRNA / Total * 100,
    tRNA = tRNA / Total * 100,
    Artifact = Artifact / Total * 100,
    Unknown = Unknown / Total * 100
  )

# Convert data to long format
mirtrace_data_primaryHBECs_percent_long <- mirtrace_data_primaryHBECs_percent %>%
  pivot_longer(cols = c(-Sample, -Total),
               names_to = 'RNA_type',
               values_to = 'Percentage') 

# Set order of samples and RNA types in the plot
mirtrace_data_primaryHBECs_percent_long <- mirtrace_data_primaryHBECs_percent_long %>%
  mutate(RNA_type = factor(RNA_type, levels = c('Unknown', 'Artifact', 'tRNA', 'rRNA', 'miRNA')), 
         Sample = factor(Sample, levels = c(primaryHBEC_untreated, 
                                            primaryHBEC_dose_0_005,
                                            primaryHBEC_dose_0_08,
                                            primaryHBEC_dose_2,
                                            primaryHBEC_dose_24)))

# Convert sample sequence ID to human-readable format
mirtrace_data_primaryHBECs_percent_long <- mirtrace_data_primaryHBECs_percent_long %>%
  mutate(
    Sample_Name = recode(
      Sample,
      "NGS-110-028_S28_R1_001" = "Untreated (1)",
      "NGS-110-016_S16_R1_001" = "Untreated (2)",
      "NGS-110-007_S7_R1_001" = "Untreated (3)",
      
      "NGS-110-029_S29_R1_001" = "Dose 0.005 (1)",
      "NGS-110-017_S17_R1_001" = "Dose 0.005 (2)",
      "NGS-110-011_S11_R1_001" = "Dose 0.005 (3)",
      
      "NGS-110-012_S12_R1_001" = "Dose 0.08 (1)",
      "NGS-110-025_S25_R1_001" = "Dose 0.08 (2)",
      "NGS-110-006_S6_R1_001" = "Dose 0.08 (3)",
      
      "NGS-110-004_S4_R1_001" = "Dose 2 (1)",
      "NGS-110-026_S26_R1_001" = "Dose 2 (2)",
      "NGS-110-003_S3_R1_001" = "Dose 2 (3)",
      
      "NGS-110-030_S30_R1_001" = "Dose 24 (1)",
      "NGS-110-019_S19_R1_001" = "Dose 24 (2)",
      "NGS-110-008_S8_R1_001" = "Dose 24 (3)",
      ))



# Create bar chart
primaryHBEC_rna_composition_plot <- ggplot(mirtrace_data_primaryHBECs_percent_long, aes(x = Percentage, y = Sample_Name, 
                                                    fill = RNA_type)) +
         geom_bar(stat = 'identity') +
  labs(x = "Proportion of Library (%)",
       y = "Sample",
       fill = 'RNA Species') +
  theme_pubclean() +
  labs_pubr() +
  theme(axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  scale_fill_manual(
    values = c(
      "miRNA"    = "steelblue4",  
      "rRNA"     = "orange2",  
      "tRNA"     = "seagreen4",  
      "Artifact" = "palevioletred3", 
      "Unknown"  = "grey85"   
    )) + 
  guides(fill = guide_legend(reverse = TRUE))


# Visualise plot
print(primaryHBEC_rna_composition_plot)


#==========================#
# Panel E - RNA Composition Plot (Calu-3)
#==========================#

# Calculate library size to enable proportion calculations 
mirtrace_data_Calu3 <- mirtrace_data_Calu3 %>%
  mutate(Total = rowSums(mirtrace_data_Calu3[-1]))

# Calculate proportion of each RNA type as percentages
mirtrace_data_calu3_percent <- mirtrace_data_Calu3 %>%
  mutate(
    miRNA = miRNA / Total * 100,
    rRNA = rRNA / Total * 100,
    tRNA = tRNA / Total * 100,
    Artifact = Artifact / Total * 100,
    Unknown = Unknown / Total * 100
  )

# Convert data to long format
mirtrace_data_calu3_percent_long <- mirtrace_data_calu3_percent %>%
  pivot_longer(cols = c(-Sample, -Total),
               names_to = 'RNA_type',
               values_to = 'Percentage') 

# Set order of samples and RNA types in the plot
mirtrace_data_calu3_percent_long <- mirtrace_data_calu3_percent_long %>%
  mutate(RNA_type = factor(RNA_type, levels = c('Unknown', 'Artifact', 'tRNA', 'rRNA', 'miRNA')), 
         Sample = factor(Sample, levels = c(calu3_untreated, 
                                            calu3_hdm_active,
                                            calu3_hdm_inactive,
                                            calu3_dog
                                            )))

# Convert sample sequence ID to human-readable format
mirtrace_data_calu3_percent_long <- mirtrace_data_calu3_percent_long %>%
  mutate(
    Sample_Name = recode(
      Sample,
      "NGS-110-021_S21_R1_001" = "Untreated (1)",
      "NGS-110-018_S18_R1_001" = "Untreated (2)",
      "NGS-110-002_S2_R1_001" = "Untreated (3)",
      
      "NGS-110-027_S27_R1_001" = "Der p 1 - Active (1)",
      "NGS-110-005_S5_R1_001" = "Der p 1 - Active (2)",
      "NGS-110-023_S23_R1_001" = "Der p 1 - Active (3)",
      
      "NGS-110-020_S20_R1_001" = "Der p 1 - Inactive (1)",
      "NGS-110-009_S9_R1_001" = "Der p 1 - Inactive (2)",
      "NGS-110-001_S1_R1_001" = "Der p 1 - Inactive (3)",
      
      "NGS-110-024_S24_R1_001" = "Can f 1 (1)",
      "NGS-110-015_S15_R1_001" = "Can f 1 (2)",
      "NGS-110-013_S13_R1_001" = "Can f 1 (3)",
    ))



# Create bar chart
calu3_rna_composition_plot <- ggplot(mirtrace_data_calu3_percent_long, aes(x = Percentage, y = Sample_Name, 
                                                    fill = RNA_type)) +
  geom_bar(stat = 'identity') +
  labs(x = "Proportion of Library (%)",
       y = "Sample",
       fill = 'RNA Species') +
  theme_pubclean() +
  labs_pubr() +
  theme(axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  scale_fill_manual(
    values = c(
      "miRNA"    = "steelblue4",  
      "rRNA"     = "orange2",  
      "tRNA"     = "seagreen4",  
      "Artifact" = "palevioletred3", 
      "Unknown"  = "grey85"   
    )) + 
  guides(fill = guide_legend(reverse = TRUE))


# Visualise plot
print(calu3_rna_composition_plot)


#==========================#
# Combine Plot to Create Figure
#==========================#

# Combine plots
figure_3 <- ((overall_mirna_boxplot | primaryHBEC_mirna_boxplot | calu3_mirna_boxplot) 
             / primaryHBEC_rna_composition_plot / calu3_rna_composition_plot) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))


# Visualise figure
print(figure_3)


# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/figure_3.pdf", plot = figure_3, height=9, width=14, dpi=600, units="in", scale=2)
