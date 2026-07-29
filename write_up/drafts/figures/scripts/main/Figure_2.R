# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(tidyr)
library(scales)
library(patchwork)


#==========================#
# Define Experimental Condition Groups
#==========================#

# Untreated
untreated <- c(
  "NGS-110-028_S28_R1_001",
  "NGS-110-016_S16_R1_001",
  "NGS-110-007_S7_R1_001"
)

# Dose 0.005
dose_0_005 <- c(
  "NGS-110-029_S29_R1_001",
  "NGS-110-017_S17_R1_001",
  "NGS-110-011_S11_R1_001"
)

# Dose 0.08
dose_0_08 <- c(
  "NGS-110-012_S12_R1_001",
  "NGS-110-025_S25_R1_001",
  "NGS-110-006_S6_R1_001"
)

# Dose 2
dose_2 <- c(
  "NGS-110-004_S4_R1_001",
  "NGS-110-026_S26_R1_001",
  "NGS-110-003_S3_R1_001"
)

# Dose 24
dose_24 <- c(
  "NGS-110-030_S30_R1_001",
  "NGS-110-019_S19_R1_001",
  "NGS-110-008_S8_R1_001"
)


#==========================#
# Panel A - Per-sequence GC Content Plot
#==========================#

# Load MultiQC data file for FastQC per-sequence GC content plot
multiqc_data_gc_content <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/FastQC_per_sequence_gc_content_plot_Percentages.txt")

# Remove the surplus text from each cell 
multiqc_data_gc_content[-1] <- lapply(multiqc_data_gc_content[-1], function(x) {
  as.numeric(sub(".*?, ", "", sub("\\)", "", x)))
})


# Convert data to long format
multiqc_data_gc_content_long <- multiqc_data_gc_content %>%
  pivot_longer(cols = -Sample,
               names_to = 'GC',
               values_to = 'Percentage') %>%
  mutate(
    GC = as.numeric(GC),
    Percentage = as.numeric(Percentage)
  )


# Add experimental condition group column
multiqc_data_gc_content_long <- multiqc_data_gc_content_long %>%
  mutate(
    Condition = case_when(
      Sample %in% untreated ~ "Untreated",
      Sample %in% dose_0_005 ~ "Dose 0.005",
      Sample %in% dose_0_08 ~ "Dose 0.08",
      Sample %in% dose_2 ~ "Dose 2",
      Sample %in% dose_24 ~ "Dose 24"))


# Create plot
GC_plot <- ggplot(multiqc_data_gc_content_long, aes(x = GC, y = Percentage, group = Sample, colour = Condition, alpha = Condition)) +
  geom_line(linewidth = 1) +
  labs(x = "GC Content (%)",
       y = "Proportion of Library (%)") +
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
  scale_alpha_manual(    
    values = c("Untreated" = 0.6,
               "Dose 0.005" = 0.6,
               "Dose 0.08" = 0.6,
               'Dose 2' = 0.6,
               'Dose 24' = 1)) +
  scale_colour_manual(
    values = c(
      'Untreated'  = 'grey60',
      'Dose 0.005' = 'skyblue2',
      'Dose 0.08'  = 'springgreen4',
      'Dose 2'     = 'orchid3',
      'Dose 24'    = 'darkorange3'
    )
  ) + 
  guides(
    colour = guide_legend(
      override.aes = list(linewidth = 5))) +
  coord_cartesian(ylim = c(0, 5)) 

# Visualise plot
print(GC_plot)



#==========================#
# Panel B - Read Length Distribution Plot
#==========================#

# Load MultiQC data file for mirtrace read length distribution plot
multiqc_read_length_data <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_length_plot.txt")

# Remove the surplus text from each cell 
multiqc_read_length_data[-1] <- lapply(multiqc_read_length_data[-1], function(x) {
  as.numeric(sub(".*?, ", "", sub("\\)", "", x)))
})

# Convert data to long format
multiqc_read_length_data_long <- multiqc_read_length_data %>%
  pivot_longer(cols = -Sample,
               names_to = 'Length',
               values_to = 'Read_Count') %>%
  mutate(
    Length = as.numeric(Length),
    Read_Count = as.numeric(Read_Count)
  )

# Add experimental condition group column
multiqc_read_length_data_long <- multiqc_read_length_data_long %>%
  mutate(
    Condition = case_when(
      Sample %in% untreated ~ "Untreated",
      Sample %in% dose_0_005 ~ "Dose 0.005",
      Sample %in% dose_0_08 ~ "Dose 0.08",
      Sample %in% dose_2 ~ "Dose 2",
      Sample %in% dose_24 ~ "Dose 24"))

# Create plot
read_length_plot <- ggplot(multiqc_read_length_data_long, aes(x = Length, y = Read_Count, group = Sample, colour = Condition, alpha = Condition)) +
  geom_line(linewidth = 1) +
  labs(x = "Read Length (bp)",
       y = "Read Count") +
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
  scale_alpha_manual(    
    values = c("Untreated" = 0.6,
               "Dose 0.005" = 0.6,
               "Dose 0.08" = 0.6,
               'Dose 2' = 0.6,
               'Dose 24' = 1)) +
  scale_colour_manual(
    values = c(
      'Untreated'  = 'grey60',
      'Dose 0.005' = 'skyblue2',
      'Dose 0.08'  = 'springgreen4',
      'Dose 2'     = 'orchid3',
      'Dose 24'    = 'darkorange3'
    )
  ) + 
  guides(
    colour = guide_legend(
      override.aes = list(linewidth = 5))) +
  coord_cartesian(ylim = c(0, 800000)) +
  scale_y_continuous(
    breaks = seq(0, 800000, by = 100000),
    labels = unit_format(scale = 1e-3, suffix = "K")
  ) 


# Visualise plot
print(read_length_plot)




#==========================#
# Panel C - Sequence Duplication Plot
#==========================#

# Load MultiQC data file for FastQC sequence duplication plot
multiqc_sequence_duplication_data <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/FastQC_sequence_duplication_levels_plot.txt")

# Rename column names to sequence duplication category
colnames(multiqc_sequence_duplication_data) <- c("Sample",
                                                 "1.0", "2.0", "3.0", "4.0", 
                                                 "5.0", "6.0", "7.0", "8.0", 
                                                 "9.0", ">10", ">50", ">100", 
                                                 ">500", ">1k", ">5k", ">10k+")


# Remove the surplus text from each cell 
multiqc_sequence_duplication_data[-1] <- lapply(multiqc_sequence_duplication_data[-1], function(x) {
  as.numeric(sub(".*?, ", "", sub("\\)", "", x)))
})


# Convert data to long format
multiqc_sequence_duplication_data_long <- multiqc_sequence_duplication_data %>%
  pivot_longer(cols = -Sample,
               names_to = 'Duplication_level',
               values_to = 'Percentage') %>%
  mutate(
    Duplication_level = Duplication_level,
    Percentage = as.numeric(Percentage)
  )


# Add experimental condition group column
multiqc_sequence_duplication_data_long <- multiqc_sequence_duplication_data_long %>%
  mutate(
    Condition = case_when(
      Sample %in% untreated ~ "Untreated",
      Sample %in% dose_0_005 ~ "Dose 0.005",
      Sample %in% dose_0_08 ~ "Dose 0.08",
      Sample %in% dose_2 ~ "Dose 2",
      Sample %in% dose_24 ~ "Dose 24"))

# Create plot
duplication_plot <- ggplot(multiqc_sequence_duplication_data_long, aes(x = factor(Duplication_level, level = colnames(multiqc_sequence_duplication_data)), 
                                                                       y = Percentage, group = Sample, colour = Condition, alpha = Condition)) +
  geom_line(linewidth = 1) +
  labs(x = "Sequence Duplication Level",
       y = "Proportion of Library (%)") +
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
  scale_alpha_manual(    
    values = c("Untreated" = 0.6,
               "Dose 0.005" = 0.6,
               "Dose 0.08" = 0.6,
               'Dose 2' = 0.6,
               'Dose 24' = 1)) +
  scale_colour_manual(
    values = c(
      'Untreated'  = 'grey60',
      'Dose 0.005' = 'skyblue2',
      'Dose 0.08'  = 'springgreen4',
      'Dose 2'     = 'orchid3',
      'Dose 24'    = 'darkorange3'
    )
  ) + 
  guides(
    colour = guide_legend(
      override.aes = list(linewidth = 5))) +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(
    breaks = seq(0, 60, by = 20)
  ) 


# Visualise plot
print(duplication_plot)



#==========================#
# Combine Plot to Create Figure
#==========================#

# Combine plots
figure_2 <- (GC_plot / read_length_plot / duplication_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))

# Visualise figure
print(figure_2)

# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_2/figure_2.pdf", plot = figure_2, height=7, width=6, dpi=600, units="in", scale=2)


