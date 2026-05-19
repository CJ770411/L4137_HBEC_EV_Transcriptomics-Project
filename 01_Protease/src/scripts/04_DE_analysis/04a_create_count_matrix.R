# Load necessary libraries
library(dplyr)
library(readr)

# Directory containing per-sample read count files
indir_all_counts <- "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/03_alignment/03e_read_quantification"

# Get paths to per-sample read count CSV files
infile_path_list <- list.files(path = indir_all_counts, pattern = "\\_trimmed_count_matrix.csv$", full.names = TRUE)

# Initialise empty list to store read counts for each sample
read_counts <- list()

# Loop through each read count file
for (file in infile_path_list) {
  
  # Read count file
  data <- read_delim(file, delim = "\t", trim_ws = TRUE)
  
  # Extract sample name from file name
  sample_name <- gsub(
    "_trimmed_count_matrix.csv", 
    "", 
    (basename(file)))
  
  # Store read counts in list
  read_counts[[sample_name]] <- data[[2]]
  
}

# Combine miRNA IDs and sample read counts into one data frame
count_matrix <- bind_cols(
  miRNA = data[[1]],
  read_counts
)

# Sum counts for mature miRNAs derived from multiple precursors
count_matrix <- aggregate(count_matrix[,2:ncol(count_matrix)], by=list(count_matrix$miRNA), FUN=sum) 


###==== edgeR Format ====###

# Create edgeR-specific count matrix
edgeR_count_matrix <- count_matrix

# Set miRNA column name as 'Symbol' for edgeR compatibility
colnames(edgeR_count_matrix)[1] <- "Symbol" 

# Save the count matrix as tab-separated TXT file
write_tsv(edgeR_count_matrix, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/edgeR_count_matrix.txt")



###==== DESeq2 Format ====###

# Create DESeq2-specific count matrix
deseq2_count_matrix <- count_matrix

# Set miRNA ID as row names
rownames(deseq2_count_matrix) <- deseq2_count_matrix[,1]

# Remove the miRNA column
deseq2_count_matrix[,1] <- NULL

# Save the count matrix as CSV file
write.csv(deseq2_count_matrix, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/deseq2_count_matrix.csv", row.names=TRUE)






