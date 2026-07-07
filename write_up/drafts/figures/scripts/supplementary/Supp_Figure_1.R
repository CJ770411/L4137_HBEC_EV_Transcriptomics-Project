library(dplyr)

# Define experimental groups
dose24_samples <- c(
  "NGS-110-008_S8_R1_001",
  "NGS-110-019_S19_R1_001",
  "NGS-110-030_S30_R1_001"
)

all_other_samples <- c(
  "NGS-110-003_S3_R1_001",
  "NGS-110-004_S4_R1_001",
  "NGS-110-006_S6_R1_001",
  "NGS-110-007_S7_R1_001",
  "NGS-110-011_S11_R1_001",
  "NGS-110-012_S12_R1_001",
  "NGS-110-016_S16_R1_001",
  "NGS-110-017_S17_R1_001",
  "NGS-110-025_S25_R1_001",
  "NGS-110-026_S26_R1_001",
  "NGS-110-028_S28_R1_001",
  "NGS-110-029_S29_R1_001"
)

#==========================#
# Read Length Distribution: mean increase
#==========================#

# Load data from 'miRTrace: Read Length Distribution'
readlength_data_protease <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/multiqc_mirtrace_length.txt")

# View data
head(readlength_data_protease)

# Rename column name to clarify it referes to length
colnames(length18_data)[colnames(length18_data) == '18'] <- 'l18'

# Isolate counts from 18bp reads
length18_data <- readlength_data_protease %>% 
  select(Sample, 'l18')

# Calculate mean counts of 18bp reads
dose24_l18_mean <- length18_data %>%
  filter(Sample %in% dose24_samples) %>%
  summarise(round(mean(l18)))

other_samples_l18_mean <- length18_data %>%
  filter(Sample %in% all_other_samples) %>%
  summarise(round(mean(l18)))

# Calculate percentage increase
percent_increase <- round((dose24_l18_mean / other_samples_l18_mean) * 100, 1)

# Output percentage increase
paste('Percentage increase for dose24 protease samples vs all other samples: ', percent_increase, '%', sep = '')



#==========================#
# Sequence Duplication: mean increase
#==========================#

# Load data from 'FastQC: Sequence Duplication Levels'
duplication_data_protease <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/fastqc_sequence_duplication_levels_plot.txt")

head(duplication_data_protease)

colnames(duplication_data_protease)[colnames(duplication_data_protease) == '15'] <- 'd15'

# Isolate counts from reads >10k duplicated
dup_10k <- duplication_data_protease %>% 
  select(Sample, 'd15')

# Clean duplication value
dup_10k <- dup_10k %>%
  mutate(d15_clean = as.numeric(str_extract(d15, "\\d+\\.\\d+"))) # takes the value containing a decimal point

# Calculate mean counts of >10k+ duplication rate
dose24_d15_mean <- dup_10k %>%
  filter(Sample %in% dose24_samples) %>%
  summarise(round(mean(d15_clean)))

other_samples_d15_mean <- dup_10k %>%
  filter(Sample %in% all_other_samples) %>%
  summarise(round(mean(d15_clean)))

# Output mean percentage duplication rates
paste("Mean percentage duplication rate for sequences duplicated >10k times in dose 24 samples: ", dose24_d15_mean, "%", sep = '')

paste("Mean percentage duplication rate for sequences duplicated >10k times in all other samples: ", other_samples_d15_mean, "%", sep = '')




#==========================#
# RNA Categories: RNA content
#==========================#

# Load data from 'miRTrace: Read Length Distribution'
rnacat_data_protease <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/02_reads_qc/trimmed/multiqc/multiqc_data/mirtrace_rna_categories_plot.txt")

rnacat_data_protease_percentages <- rnacat_data_protease %>%
  mutate(total = rowSums(rnacat_data_protease[-1]),
         rRNA_percent = round((rRNA / total) * 100, 2),
         tRNA_percent = round((tRNA / total) * 100, 2),
         miRNA_percent = round((miRNA / total) * 100, 2))

# Calculate mean counts of >10k+ duplication rate
dose24_rna_content_mean <- rnacat_data_protease_percentages %>%
  filter(Sample %in% dose24_samples) %>%
  summarise(mean(miRNA_percent),
            mean(rRNA_percent),
            mean(tRNA_percent))

other_samples_rna_content_mean <- rnacat_data_protease_percentages %>%
  filter(Sample %in% all_other_samples) %>%
  summarise(mean(miRNA_percent),
            mean(rRNA_percent),
            mean(tRNA_percent))

# Output mean percentage RNA contents
print("Mean percentage RNA content in dose 24 protease samples:" )
print(dose24_rna_content_mean)

print("Mean percentage RNA content in all other samples protease samples:" )
print(other_samples_rna_content_mean)
