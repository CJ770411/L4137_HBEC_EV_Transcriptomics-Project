library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

# Load data from miRTrace 'RNA Categories' plot in multiQC
mirtrace_data_protease <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_rna_categories_plot.txt")
mirtrace_data_hdm_dog <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_rna_categories_plot.txt")

# Define experimental groups (Protease)
primary_protease_untreated <- c(
  "NGS-110-028_S28_R1_001",
  "NGS-110-016_S16_R1_001",
  "NGS-110-007_S7_R1_001"
)

primary_dose_0005 <- c(
  "NGS-110-029_S29_R1_001",
  "NGS-110-017_S17_R1_001",
  "NGS-110-011_S11_R1_001"
)

primary_dose_008 <- c(
  "NGS-110-012_S12_R1_001",
  "NGS-110-025_S25_R1_001",
  "NGS-110-006_S6_R1_001"
)

primary_dose_2 <- c(
  "NGS-110-004_S4_R1_001",
  "NGS-110-026_S26_R1_001",
  "NGS-110-003_S3_R1_001"
)

primary_dose_24 <- c(
  "NGS-110-030_S30_R1_001",
  "NGS-110-019_S19_R1_001",
  "NGS-110-008_S8_R1_001"
)

# Define experimental groups (HDM_Dog)
calu3_hdm_active <- c(
  "NGS-110-027_S27_R1_001",
  "NGS-110-005_S5_R1_001",
  "NGS-110-023_S23_R1_001"
)

calu3_hdm_inactive <- c(
  "NGS-110-020_S20_R1_001",
  "NGS-110-009_S9_R1_001",
  "NGS-110-001_S1_R1_001"
)

calu3_dog <- c(
  "NGS-110-024_S24_R1_001",
  "NGS-110-015_S15_R1_001",
  "NGS-110-013_S13_R1_001"
)

calu3_poly_IC <- c(
  "NGS-110-014_S14_R1_001",
  "NGS-110-010_S10_R1_001",
  "NGS-110-022_S22_R1_001"
)

calu3_control <- c(
  "NGS-110-021_S21_R1_001",
  "NGS-110-018_S18_R1_001",
  "NGS-110-002_S2_R1_001"
)

#==========================#
# Calculating mean miRNA content (treated vs untreated)
#==========================#


# Calculate mean values (Protease)
mean_mirna_content_protease_untreated <- mirtrace_data_protease %>%
  filter(Sample %in% primary_protease_untreated) %>%
  summarise(round(mean(miRNA)))

mean_mirna_content_protease_treated <- mirtrace_data_protease %>%
  filter(Sample %in% c(primary_dose_0005,
                       primary_dose_008,
                       primary_dose_2,
                       primary_dose_24)) %>%
  summarise(round(mean(miRNA)))



# Calculate mean values (HDM_Dog)
mean_mirna_content_hdm_dog_untreated <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_control) %>%
  summarise(round(mean(miRNA)))

mean_mirna_content_hdm_dog_treated <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% c(calu3_hdm_active,
                       calu3_hdm_inactive,
                       calu3_dog,
                       calu3_poly_IC)) %>%
  summarise(round(mean(miRNA)))



# Print mean values
print(paste("Mean miRNA content in HDM_Dog - UNTREATED:",
            mean_mirna_content_hdm_dog_untreated))
print(paste("Mean miRNA content in HDM_Dog - TREATED:",
            mean_mirna_content_hdm_dog_treated))

print(paste("Mean miRNA content in Protease - UNTREATED:",
            mean_mirna_content_protease_untreated))
print(paste("Mean miRNA content in Protease - TREATED:",
            mean_mirna_content_protease_treated))



#==========================#
# Boxplot
#==========================#

# Extract miRNA counts and create condition column for all contrasts

# Create condition columns (Protease)

# All protease samples
mirna_prim_protease <- mirtrace_data_protease %>%
  mutate(Condition = 'Protease') %>%
  select(miRNA, Sample, Condition) 

# Conditions
mirna_primary_protease_0005 <- mirtrace_data_protease %>%
  filter(Sample %in% primary_dose_0005) %>%
  mutate(Condition = 'Dose 0.005') %>%
  select(miRNA, Sample, Condition)

mirna_primary_protease_008 <- mirtrace_data_protease %>%
  filter(Sample %in% primary_dose_008) %>%
  mutate(Condition = 'Dose 0.08') %>%
  select(miRNA, Sample, Condition)

mirna_primary_protease_2 <- mirtrace_data_protease %>%
  filter(Sample %in% primary_dose_2) %>%
  mutate(Condition = 'Dose 2') %>%
  select(miRNA, Sample, Condition)

mirna_primary_protease_24 <- mirtrace_data_protease %>%
  filter(Sample %in% primary_dose_24) %>%
  mutate(Condition = 'Dose 24') %>%
  select(miRNA, Sample, Condition)

mirna_primary_protease_untreated <- mirtrace_data_protease %>%
  filter(Sample %in% primary_protease_untreated) %>%
  mutate(Condition = 'Untreated') %>%
  select(miRNA, Sample, Condition)


# Create condition columns (HDM_Dog)

# All HDM_Dog samples
mirna_calu3_hdm_dog <- mirtrace_data_hdm_dog %>%
  mutate(Condition = 'HDM_Dog') %>%
  select(miRNA, Sample, Condition) 

# Conditions
mirna_calu3_hdm_active <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_hdm_active) %>%
    mutate(Condition = 'HDM_Active') %>%
    select(miRNA, Sample, Condition)

mirna_calu3_hdm_inactive <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_hdm_inactive) %>%
  mutate(Condition = 'HDM_Inactive') %>%
  select(miRNA, Sample, Condition)

mirna_calu3_dog <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_dog) %>%
  mutate(Condition = 'Dog') %>%
  select(miRNA, Sample, Condition)

mirna_calu3_Poly_IC <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_poly_IC) %>%
  mutate(Condition = 'Poly_IC') %>%
  select(miRNA, Sample, Condition)

mirna_calu3_untreated <- mirtrace_data_hdm_dog %>%
  filter(Sample %in% calu3_control) %>%
  mutate(Condition = 'Untreated') %>%
  select(miRNA, Sample, Condition)




# Combine data into one dataframe
boxplot_data_protease <- rbind(mirna_primary_protease_untreated,
                               mirna_primary_protease_0005,
                               mirna_primary_protease_008,
                               mirna_primary_protease_2,
                               mirna_primary_protease_24)

boxplot_data_hdm_dog <- rbind(mirna_calu3_untreated,
                              mirna_calu3_hdm_active,
                              mirna_calu3_hdm_inactive,
                              mirna_calu3_dog,
                              mirna_calu3_Poly_IC)


boxplot_data_all <- rbind(mirna_prim_protease,
                              mirna_calu3_hdm_dog)


# Plot boxplot
ggplot(boxplot_data_protease, aes(x = Condition, y = miRNA, fill = Condition, alpha = 0.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  theme_pubr(legend = 'right') +
  labs_pubr() 

ggplot(boxplot_data_hdm_dog, aes(x = Condition, y = miRNA, fill = Condition, alpha = 0.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  theme_pubr(legend = 'right') +
  labs_pubr() 

ggplot(boxplot_data_all, aes(x = Condition, y = miRNA, fill = Condition, alpha = 0.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 8,
               size = 2, color = "red") +
  scale_fill_viridis_d() +
  theme_pubr(legend = 'right') +
  labs_pubr() 

#[need to decide what the boxplot should actually show, e.g., untreated vs treated, contrasts or only two boxes with primary cells vs calu-3]


#==========================#
# Library Size
#==========================#

# Load data from 'MultiQC: General Statistics'
generalstats_data_protease <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/multiqc_general_stats.txt")
generalstats_data_hdm_dog <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/02_reads_qc/trimmed/multiqc/multiqc_data/multiqc_general_stats.txt")

# Calculate mean library size
mean_library_size_protease <- generalstats_data_protease %>%
  summarise(mean(`fastqc-total_sequences`) * 1000000)

mean_library_size_hdm_dog <- generalstats_data_hdm_dog %>%
  summarise(mean(`fastqc-total_sequences`) * 1000000)

# Output mmean library size 
print("Mean library size for protease samples:" )
print(mean_library_size_protease)

print("Mean library size for HDM_Dog samples:" )
print(mean_library_size_hdm_dog)






