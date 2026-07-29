# Load necessary libraries
library(dplyr)
library(readr)

# Directory containing per-sample read count files
indir_all_counts <- "../../../results/primaryHBEC/05_read_quantification"

# Get paths to per-sample read count CSV files
infile_path_list <- list.files(path = indir_all_counts, pattern = "\\_trimmed_count_matrix.csv$", full.names = TRUE)

# Initialise empty list to store read counts for each sample
read_counts <- list()
# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix
count_matrix <- read_tsv("../../../results/primaryHBEC/06_create_count_matrix/edgeR_count_matrix.txt")

# Create experimental condition labels 
group <- factor(c("Control", "Control", "Control",
                  "Dose_0_005", "Dose_0_005", "Dose_0_005",
                  "Dose_0_08", "Dose_0_08", "Dose_0_08",
                  "Dose_2", "Dose_2", "Dose_2",
                  "Dose_24", "Dose_24", "Dose_24"))

# Set order of experimental conditions
group <- factor(group,
                levels=c("Control",
                         "Dose_0_005",
                         "Dose_0_08",
                         "Dose_2",
                         "Dose_24"))

# Create DGEList object containing read data and condition
dge <- DGEList(counts=count_matrix, group=group)

# Sanity check of DGEList structure
dge$samples

# Observe number of miRNAs pre-filtering
nrow(dge)

# Identify miRNA with sufficient counts in samples, n, where n=sample count in smallest group
keep <- filterByExpr(dge)

# Filter miRNA to keep only miRNAs with worthwhile counts
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Observe number of miRNAs post-filtering
nrow(dge) 

# Perform TMM normalisation
dge <- normLibSizes(dge)

# Visualise clustering of samples to detect outliers
plotMDS(dge, labels = dge$samples$group, col=rep(c("red", "blue", "darkorange", "darkgreen", "purple"), each = 3))

# Create design matrix defining pairwise comparisons
design <- model.matrix(~group)

# Observe design matrix
design

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Visualise biological coefficient of variation (BCV) 
plotBCV(dge)

# Fit quasi-likelihood F test model for each miRNA
fit <- glmQLFit(dge, design, robust = TRUE)

# Save 'dge' and 'group' objects for PCA plot generation
saveRDS(
  dge,
  "../../../results/primaryHBEC/07_DE_analysis/dge.rds"
)

saveRDS(
  group,
  "../../../results/primaryHBEC/07_DE_analysis/dge_group.rds"
)



#==========================#
# Differential Expression      
#==========================#


logFC_threshold <- 1 # miRNA expression threshold

###==== Dose: 0.005 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose0005 <- glmQLFTest(fit, coef="groupDose_0_005")

# Calculate differential expression
all_results_dose0005 <- as.data.frame(topTags(qlf_dose0005, n=nrow(qlf_dose0005)))

# Extract significant results 
sig_results_FDR005_dose0005 <- all_results_dose0005[all_results_dose0005$FDR<0.05 & abs(all_results_dose0005$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_dose0005, "../../../results/primaryHBEC/07_DE_analysis/all_results_dose0005.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose0005, "../../../results/primaryHBEC/07_DE_analysis/sig_results_FDR005_dose0005.csv", row.names = FALSE)



###==== Dose: 0.08 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose008 <- glmQLFTest(fit, coef="groupDose_0_08")

# Calculate differential expression
all_results_dose008 <- as.data.frame(topTags(qlf_dose008, n=nrow(qlf_dose008)))

# Extract significant results 
sig_results_FDR005_dose008 <- all_results_dose008[all_results_dose008$FDR<0.05 & abs(all_results_dose008$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_dose008, "../../../results/primaryHBEC/07_DE_analysis/all_results_dose008.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose008, "../../../results/primaryHBEC/07_DE_analysis/sig_results_FDR005_dose008.csv", row.names = FALSE)



###==== Dose: 2 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose2 <- glmQLFTest(fit, coef="groupDose_2")

# Calculate differential expression
all_results_dose2 <- as.data.frame(topTags(qlf_dose2, n=nrow(qlf_dose2)))

# Extract significant results 
sig_results_FDR005_dose2 <- all_results_dose2[all_results_dose2$FDR<0.05 & abs(all_results_dose2$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_dose2, "../../../results/primaryHBEC/07_DE_analysis/all_results_dose2.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose2, "../../../results/primaryHBEC/07_DE_analysis/sig_results_FDR005_dose2.csv", row.names = FALSE)



###==== Dose: 24 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose24 <- glmQLFTest(fit, coef="groupDose_24")

# Calculate differential expression
all_results_dose24 <- as.data.frame(topTags(qlf_dose24, n=nrow(qlf_dose24)))

# Extract significant results 
sig_results_FDR005_dose24 <- all_results_dose24[all_results_dose24$FDR<0.05 & abs(all_results_dose24$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_dose24, "../../../results/primaryHBEC/07_DE_analysis/all_results_dose24.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose24, "../../../results/primaryHBEC/07_DE_analysis/sig_results_FDR005_dose24.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_dose24,
  "../../../results/primaryHBEC/07_DE_analysis/all_results_dose24.rds"
)



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
write_tsv(edgeR_count_matrix, "../../../results/primaryHBEC/06_create_count_matrix/edgeR_count_matrix.txt")



###==== DESeq2 Format ====###

# Create DESeq2-specific count matrix
deseq2_count_matrix <-as.data.frame(count_matrix)

# Set miRNA ID as row names
rownames(deseq2_count_matrix) <- deseq2_count_matrix[,1]

# Remove the miRNA column
deseq2_count_matrix[,1] <- NULL

# Save the count matrix as CSV file
write.csv(deseq2_count_matrix, "../../../results/primaryHBEC/06_create_count_matrix/deseq2_count_matrix.csv", row.names=TRUE)






