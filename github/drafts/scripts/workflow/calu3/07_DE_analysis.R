# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix
count_matrix <- read_tsv("../../../results/calu3/06_create_count_matrix/edgeR_count_matrix.txt")

# Create experimental condition labels 
group <- factor(c("Der_p1_Active", "Der_p1_Active", "Der_p1_Active",
                  "Can_F_1", "Can_F_1", "Can_F_1",
                  "Der_p1_Inactive", "Der_p1_Inactive", "Der_p1_Inactive",
                  "Unstimulated", "Unstimulated", "Unstimulated"))

# Set order of experimental conditions
group <- factor(
  group,
  levels = c(
    "Unstimulated",
    "Der_p1_Active",
    "Can_F_1",
    "Der_p1_Inactive"
  )
)

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
plotMDS(dge, labels = dge$samples$group, col=rep(c("red", "blue", "darkorange", "darkgreen"), each = 3))

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
  "../../../results/calu3/07_DE_analysis/dge.rds"
)

saveRDS(
  group,
  "../../../results/calu3/07_DE_analysis/dge_group.rds"
)



#==========================#
# Differential Expression      
#==========================#

logFC_threshold <- 1 # miRNA expression threshold


###==== Der_p1_Active ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Der_p1_Active <- glmQLFTest(fit, coef="groupDer_p1_Active")

# Calculate differential expression
all_results_qlf_Der_p1_Active <- as.data.frame(topTags(qlf_Der_p1_Active, n=nrow(qlf_Der_p1_Active)))

# Extract significant results 
sig_results_FDR005_qlf_Der_p1_Active <- all_results_qlf_Der_p1_Active[all_results_qlf_Der_p1_Active$FDR<0.05 & abs(all_results_qlf_Der_p1_Active$logFC) >= logFC_threshold,]
sig_results_FDR02_qlf_Der_p1_Active <- all_results_qlf_Der_p1_Active[all_results_qlf_Der_p1_Active$FDR<0.2 & abs(all_results_qlf_Der_p1_Active$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Der_p1_Active, "../../../results/calu3/07_DE_analysis/all_results_qlf_Der_p1_Active.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Der_p1_Active, "../../../results/calu3/07_DE_analysis/sig_results_FDR005_qlf_Der_p1_Active.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Der_p1_Active, "../../../results/calu3/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Active.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Der_p1_Active,
  "../../../results/calu3/07_DE_analysis/all_results_qlf_Der_p1_Active.rds"
)


###==== Can_F_1 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Can_F_1 <- glmQLFTest(fit, coef="groupCan_F_1")

# Calculate differential expression
all_results_qlf_Can_F_1 <- as.data.frame(topTags(qlf_Can_F_1, n=nrow(qlf_Can_F_1)))

# Extract significant results 
sig_results_FDR005_qlf_Can_F_1 <- all_results_qlf_Can_F_1[all_results_qlf_Can_F_1$FDR<0.05 & abs(all_results_qlf_Can_F_1$logFC) >= logFC_threshold,]
sig_results_FDR02_qlf_Can_F_1 <- all_results_qlf_Can_F_1[all_results_qlf_Can_F_1$FDR<0.2 & abs(all_results_qlf_Can_F_1$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Can_F_1, "../../../results/calu3/07_DE_analysis/all_results_qlf_Can_F_1.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Can_F_1, "../../../results/calu3/07_DE_analysis/sig_results_FDR005_qlf_Can_F_1.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Can_F_1, "../../../results/calu3//07_DE_analysis/sig_results_FDR02_qlf_Can_F_1.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Can_F_1,
  "../../../results/calu3/07_DE_analysis/all_results_qlf_Can_F_1.rds"
)


###==== Der_p1_Inactive ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Der_p1_Inactive <- glmQLFTest(fit, coef="groupDer_p1_Inactive")

# Calculate differential expression
all_results_qlf_Der_p1_Inactive <- as.data.frame(topTags(qlf_Der_p1_Inactive, n=nrow(qlf_Der_p1_Inactive)))

# Extract significant results 
sig_results_FDR005_qlf_Der_p1_Inactive <- all_results_qlf_Der_p1_Inactive[all_results_qlf_Der_p1_Inactive$FDR<0.05 & abs(all_results_qlf_Der_p1_Inactive$logFC) >= logFC_threshold,]
sig_results_FDR02_qlf_Der_p1_Inactive <- all_results_qlf_Der_p1_Inactive[all_results_qlf_Der_p1_Inactive$FDR<0.2 & abs(all_results_qlf_Der_p1_Inactive$logFC) >= logFC_threshold,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Der_p1_Inactive, "../../../results/calu3/07_DE_analysis/all_results_qlf_Der_p1_Inactive.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Der_p1_Inactive, "../../../results/calu3/07_DE_analysis/sig_results_FDR005_qlf_Der_p1_Inactive.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Der_p1_Inactive, "../../../results/calu3/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Inactive.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Der_p1_Inactive,
  "../../../results/calu3/07_DE_analysis/all_results_qlf_Der_p1_Inactive.rds"
)


