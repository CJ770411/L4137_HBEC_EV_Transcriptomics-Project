# Load necessary libraries
library(dplyr)
library(readr)

# Directory containing per-sample read count files
indir_all_counts <- "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/05_read_quantification"

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
count_matrix <- count_matrix[, c("Symbol",
                                 "NGS-110-028_S28_R1_001", "NGS-110-016_S16_R1_001", "NGS-110-007_S7_R1_001",
                                 "NGS-110-029_S29_R1_001", "NGS-110-017_S17_R1_001", "NGS-110-011_S11_R1_001", 
                                 "NGS-110-012_S12_R1_001", "NGS-110-025_S25_R1_001", "NGS-110-006_S6_R1_001", 
                                 "NGS-110-004_S4_R1_001", "NGS-110-026_S26_R1_001", "NGS-110-003_S3_R1_001", 
                                 "NGS-110-030_S30_R1_001", "NGS-110-019_S19_R1_001", "NGS-110-008_S8_R1_001")]

# Set new column names
colnames(count_matrix) <- c("Symbol", 
                            "Control__1", "Control__2", "Control__3",
                            "Dose0_005__1", "Dose0_005__2", "Dose0_005__3",
                            "Dose0_08__1", "Dose0_08__2", "Dose0_08__3",
                            "Dose2__1", "Dose2__2", "Dose2__3",
                            "Dose24__1", "Dose24__2", "Dose24__3")

# Remove rows with 0 counts
count_matrix <- count_matrix[rowSums(count_matrix[,-1]) > 0, ] # excludes miRNA column (col1) because it doesn't work on character vectors


###==== edgeR Format ====###

# Create edgeR-specific count matrix
edgeR_count_matrix <- count_matrix

# Save the count matrix as tab-separated TXT file
write_tsv(edgeR_count_matrix, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/06_create_count_matrix/edgeR_count_matrix.txt")



###==== DESeq2 Format ====###

# Create DESeq2-specific count matrix
deseq2_count_matrix <- as.data.frame(count_matrix)

# Set miRNA ID as row names
rownames(deseq2_count_matrix) <- deseq2_count_matrix[,1]

# Remove the miRNA column
deseq2_count_matrix[,1] <- NULL

# Save the count matrix as CSV file
write.csv(deseq2_count_matrix, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_03/06_create_count_matrix/deseq2_count_matrix.csv", row.names=TRUE)






