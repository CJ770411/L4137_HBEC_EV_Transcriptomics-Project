#if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
#BiocManager::install("edgeR")

# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix from '04b_DE_analysis_edgeR.R'
count_matrix <- read_tsv("/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04a_create_count_matrix/edgeR_count_matrix.txt")

# Convert count matrix to data frame
count_matrix <- as.data.frame(count_matrix)

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
fit <- glmQLFit(dge, design)


#==========================#
# Differential Expression      
#==========================#


###==== Dose: 0.005 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose0005 <- glmQLFTest(fit, coef="groupDose_0_005")

# Calculate differential expression
all_results_dose0005 <- as.data.frame(topTags(qlf_dose0005, n=nrow(qlf_dose0005)))

# Extract significant results 
sig_results_FDR005_dose0005 <- all_results_dose0005[all_results_dose0005$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose0005, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/all_results_dose0005.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose0005, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose0005.csv", row.names = FALSE)



###==== Dose: 0.08 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose008 <- glmQLFTest(fit, coef="groupDose_0_08")

# Calculate differential expression
all_results_dose008 <- as.data.frame(topTags(qlf_dose008, n=nrow(qlf_dose008)))

# Extract significant results 
sig_results_FDR005_dose008 <- all_results_dose008[all_results_dose008$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose008, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/all_results_dose008.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose008, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose008.csv", row.names = FALSE)



###==== Dose: 2 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose2 <- glmQLFTest(fit, coef="groupDose_2")

# Calculate differential expression
all_results_dose2 <- as.data.frame(topTags(qlf_dose2, n=nrow(qlf_dose2)))

# Extract significant results 
sig_results_FDR005_dose2 <- all_results_dose2[all_results_dose2$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose2, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/all_results_dose2.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose2, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose2.csv", row.names = FALSE)



###==== Dose: 24 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose24 <- glmQLFTest(fit, coef="groupDose_24")

# Calculate differential expression
all_results_dose24 <- as.data.frame(topTags(qlf_dose24, n=nrow(qlf_dose24)))

# Extract significant results 
sig_results_FDR005_dose24 <- all_results_dose24[all_results_dose24$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose24, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/all_results_dose24.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose24, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/sig_results_FDR005_dose24.csv", row.names = FALSE)


# [test code]
plotQLDisp(fit)
barplot(dge$samples$lib.size)
logcpm <- cpm(dge, log=TRUE)
write.csv(logcpm_df, "/Users/christopherjanschke/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/methods_sections/04_DE_analysis/04b_DE_analysis_edgeR/logcpm.csv", row.names = TRUE)
colnames(dge$genes$Symbol)
class(as.list(dge$genes$Symbol))


logcpm_df <- as.data.frame(logcpm)
rownames(logcpm_df) <- as.list(dge$genes$Symbol)

results_0_005 %>% head()


