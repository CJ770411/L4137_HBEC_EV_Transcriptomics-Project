# Load necessary libraries
library(dplyr)
library(readr)

# Path to directory containing read counts for individual samples
indir_all_counts <- "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/03_alignment/03e_read_quantification"

# Create list of paths to CSV files in the 'indir_all_counts' directory
infile_path_list <- list.files(path = indir_all_counts, pattern = "\\_trimmed_count_matrix.csv$", full.names = TRUE)

# Create empty list which will become the count matrix
read_counts <- list()

# Loop over each file and add the read counts into the count matrix
for (file in infile_path_list) {
  # Read in .csv file containing read counts
  data <- read_delim(file, delim = "\t", trim_ws = TRUE)
  
  # Extract column 2 containing read counts
  count_col <- data[,2]
  
  # Extract sample name
  sample_name <- gsub("_trimmed_count_matrix.csv", "", (basename(file)))
  
  # Add sample name as column name then add read counts to count matrix
  read_counts[[sample_name]] <- data[[2]]
  
}

count_matrix <- bind_cols(read_counts)

count_matrix <- bind_cols(
  miRNA = data[[1]],
  count_matrix
)


write.csv(count_matrix, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/03_alignment/03e_read_quantification/miRDeep2_count_matrix.csv", row.names=FALSE)


rm(list=ls())



