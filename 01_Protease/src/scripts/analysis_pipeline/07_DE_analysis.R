#if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
#BiocManager::install("edgeR")

# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix from '04a_create_count_matrix'
count_matrix <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/06_create_count_matrix/edgeR_count_matrix.txt")

# Convert count matrix to data frame
#count_matrix <- as.data.frame(count_matrix)

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
fit <- glmQLFit(dge, design, robust = TRUE)

# Save 'dge' and 'group' objects for PCA plot generation
saveRDS(
  dge,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/dge.rds"
)

saveRDS(
  group,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/dge_group.rds"
)



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
write.csv(all_results_dose0005, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose0005.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose0005, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose0005.csv", row.names = FALSE)



###==== Dose: 0.08 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose008 <- glmQLFTest(fit, coef="groupDose_0_08")

# Calculate differential expression
all_results_dose008 <- as.data.frame(topTags(qlf_dose008, n=nrow(qlf_dose008)))

# Extract significant results 
sig_results_FDR005_dose008 <- all_results_dose008[all_results_dose008$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose008, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose008.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose008, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose008.csv", row.names = FALSE)



###==== Dose: 2 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose2 <- glmQLFTest(fit, coef="groupDose_2")

# Calculate differential expression
all_results_dose2 <- as.data.frame(topTags(qlf_dose2, n=nrow(qlf_dose2)))

# Extract significant results 
sig_results_FDR005_dose2 <- all_results_dose2[all_results_dose2$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose2, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose2.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose2, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose2.csv", row.names = FALSE)



###==== Dose: 24 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_dose24 <- glmQLFTest(fit, coef="groupDose_24")

# Calculate differential expression
all_results_dose24 <- as.data.frame(topTags(qlf_dose24, n=nrow(qlf_dose24)))

# Extract significant results 
sig_results_FDR005_dose24 <- all_results_dose24[all_results_dose24$FDR<0.05,]

# Save the DE results as CSV file
write.csv(all_results_dose24, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose24.csv", row.names = FALSE)
write.csv(sig_results_FDR005_dose24, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_dose24.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_dose24,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose24.rds"
)


