# Load necessary libraries
library(dplyr)
library(readr)

# Directory containing per-sample read count files
indir_all_counts <- "../../../results/calu3/05_read_quantification/all_samples"

# Get paths to per-sample read count CSV files
infile_path_list <- list.files(path = indir_all_counts, pattern = "\\_trimmed_count_matrix.csv$", full.names = TRUE)

# Initialise empty list to store read counts for each sample
read_counts <- list()

# Loop through each read count file
for (file in infile_path_list) {
  
  # Read count file
  data <- read_delim(file, delim = "\t", trim_ws = TRUE) %>%
    dplyr::group_by(`#miRNA`) %>%
    
    # Get maximum mature miRNA count if multiple precursors
    dplyr::summarise(read_count = max(read_count), .groups = "drop") 
  
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

# Set miRNA column name as 'Symbol' for clarity
colnames(count_matrix)[1] <- "Symbol" 

# Re-order columns to match condition/replicate order
count_matrix <- count_matrix[, c(
  "Symbol",
  "NGS-110-027_S27_R1_001",
  "NGS-110-005_S5_R1_001",
  "NGS-110-023_S23_R1_001",
  "NGS-110-024_S24_R1_001",
  "NGS-110-015_S15_R1_001",
  "NGS-110-013_S13_R1_001",
  "NGS-110-020_S20_R1_001",
  "NGS-110-009_S9_R1_001",
  "NGS-110-001_S1_R1_001",
  "NGS-110-021_S21_R1_001",
  "NGS-110-018_S18_R1_001",
  "NGS-110-002_S2_R1_001"
)]

# Set new column names
colnames(count_matrix) <- c("Symbol", 
                            "Der_p1_Active__1", "Der_p1_Active__2", "Der_p1_Active__3",
                            "Can_F_1__1", "Can_F_1__2", "Can_F_1__3",
                            "Der_p1_Inactive__1", "Der_p1_Inactive__2", "Der_p1_Inactive__3",
                            "Unstimulated__1", "Unstimulated__2", "Unstimulated__3")

# Remove rows with 0 counts
count_matrix <- count_matrix[rowSums(count_matrix[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors


###==== edgeR Format ====###

# Create edgeR-specific count matrix
edgeR_count_matrix <- count_matrix

# Save the count matrix as tab-separated TXT file
write_tsv(edgeR_count_matrix, "../../../results/calu3/06_create_count_matrix/edgeR_count_matrix.txt")



###==== DESeq2 Format ====###
# Note: This format was no used in the workflow.

# Create DESeq2-specific count matrix
deseq2_count_matrix <- as.data.frame(count_matrix)

# Set miRNA ID as row names
rownames(deseq2_count_matrix) <- deseq2_count_matrix[,1]

# Remove the miRNA column
deseq2_count_matrix[,1] <- NULL

# Save the count matrix as CSV file
write.csv(deseq2_count_matrix, "../../../results/calu3/06_create_count_matrix/deseq2_count_matrix.csv", row.names=TRUE)






