#if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
#BiocManager::install("edgeR")

# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix from '04a_create_count_matrix.R'
count_matrix <- read_tsv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/edgeR_count_matrix.txt")

# Convert count matrix to data frame
count_matrix <- as.data.frame(count_matrix)

# Create experimental condition labels 
group <- factor(c("Control", "Control", "Control",
                  "Dose_0.005", "Dose_0.005", "Dose_0.005",
                  "Dose_0.08", "Dose_0.08", "Dose_0.08",
                  "Dose_2", "Dose_2", "Dose_2",
                  "Dose_24", "Dose_24", "Dose_24"))

# Set order of experimental conditions
group <- factor(group,
                levels=c("Control",
                         "Dose_0.005",
                         "Dose_0.08",
                         "Dose_2",
                         "Dose_24"))

# Re-order columns to match condition label order
count_matrix <- count_matrix[, c("Symbol",
                                 "NGS-110-028_S28_R1_001", "NGS-110-016_S16_R1_001", "NGS-110-007_S7_R1_001",
                                 "NGS-110-029_S29_R1_001", "NGS-110-017_S17_R1_001", "NGS-110-011_S11_R1_001", 
                                 "NGS-110-012_S12_R1_001", "NGS-110-025_S25_R1_001", "NGS-110-006_S6_R1_001", 
                                 "NGS-110-004_S4_R1_001", "NGS-110-026_S26_R1_001", "NGS-110-003_S3_R1_001", 
                                 "NGS-110-030_S30_R1_001", "NGS-110-019_S19_R1_001", "NGS-110-008_S8_R1_001")]

# Create DGEList object containing read data and condition
dge <- DGEList(counts=count_matrix, group=group)

# Sanity check of DGEList structure
dge$samples

# Observe number of miRNAs pre-filtering
nrow(dge)

# Identify miRNA with sufficient counts in samples, n, where n=sample count in smallest group.
keep <- filterByExpr(dge)

# Filter miRNA to keep only miRNAs with worthwhile counts
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Observe number of miRNAs post-filtering
nrow(dge) 

# Perform TMM normalisation
dge <- normLibSizes(dge)

# Visualise clustering of samples to detect outliers
plotMDS(dge, labels = dge$samples[,1])

# Create design matrix defining pairwise comparisons
design <- model.matrix(~group)

# Observe design matrix
design

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Visualise biological coefficient of variation (BCV) 
plotBCV(dge)

fit <- glmQLFit(dge, design)

# Quasi-likelihood
qlf_0.005 <- glmQLFTest(fit, coef="groupDose_0.005")

results_0.005 <- as.data.frame(topTags(qlf_0.005, n=nrow(qlf_0.005)))

dim(results_0.005)

write.csv(results_0.005, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/results_0.005.csv", row.names = FALSE)

