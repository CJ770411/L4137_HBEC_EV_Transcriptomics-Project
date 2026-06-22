#if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
#BiocManager::install("edgeR")

# Load necessary libraries
library(edgeR)
library(readr)
library(dplyr)

# Read in count matrix
count_matrix <- read_tsv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/06_create_count_matrix/edgeR_count_matrix.txt")

# Create experimental condition labels 
group <- factor(c("Der_p1_Active", "Der_p1_Active", "Der_p1_Active",
                  "Can_F_1", "Can_F_1", "Can_F_1",
                  "Der_p1_Inactive", "Der_p1_Inactive", "Der_p1_Inactive",
                  "Poly_IC", "Poly_IC", "Poly_IC",
                  "Unstimulated", "Unstimulated", "Unstimulated"))

# Set order of experimental conditions
group <- factor(
  group,
  levels = c(
    "Unstimulated",
    "Der_p1_Active",
    "Can_F_1",
    "Der_p1_Inactive",
    "Poly_IC"
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
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/dge.rds"
)

saveRDS(
  group,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/dge_group.rds"
)



#==========================#
# Differential Expression      
#==========================#


###==== Der_p1_Active ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Der_p1_Active <- glmQLFTest(fit, coef="groupDer_p1_Active")

# Calculate differential expression
all_results_qlf_Der_p1_Active <- as.data.frame(topTags(qlf_Der_p1_Active, n=nrow(qlf_Der_p1_Active)))

# Extract significant results 
sig_results_FDR005_qlf_Der_p1_Active <- all_results_qlf_Der_p1_Active[all_results_qlf_Der_p1_Active$FDR<0.05,]
sig_results_FDR02_qlf_Der_p1_Active <- all_results_qlf_Der_p1_Active[all_results_qlf_Der_p1_Active$FDR<0.2,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Der_p1_Active, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Der_p1_Active.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Der_p1_Active, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_qlf_Der_p1_Active.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Der_p1_Active, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Active.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Der_p1_Active,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Der_p1_Active.rds"
)


###==== Can_F_1 ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Can_F_1 <- glmQLFTest(fit, coef="groupCan_F_1")

# Calculate differential expression
all_results_qlf_Can_F_1 <- as.data.frame(topTags(qlf_Can_F_1, n=nrow(qlf_Can_F_1)))

# Extract significant results 
sig_results_FDR005_qlf_Can_F_1 <- all_results_qlf_Can_F_1[all_results_qlf_Can_F_1$FDR<0.05,]
sig_results_FDR02_qlf_Can_F_1 <- all_results_qlf_Can_F_1[all_results_qlf_Can_F_1$FDR<0.2,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Can_F_1, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Can_F_1.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Can_F_1, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_qlf_Can_F_1.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Can_F_1, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR02_qlf_Can_F_1.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Can_F_1,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Can_F_1.rds"
)


###==== Der_p1_Inactive ====###

# Perform quasi-likelihood F test on each miRNA
qlf_Der_p1_Inactive <- glmQLFTest(fit, coef="groupDer_p1_Inactive")

# Calculate differential expression
all_results_qlf_Der_p1_Inactive <- as.data.frame(topTags(qlf_Der_p1_Inactive, n=nrow(qlf_Der_p1_Inactive)))

# Extract significant results 
sig_results_FDR005_qlf_Der_p1_Inactive <- all_results_qlf_Der_p1_Inactive[all_results_qlf_Der_p1_Inactive$FDR<0.05,]
sig_results_FDR02_qlf_Der_p1_Inactive <- all_results_qlf_Der_p1_Inactive[all_results_qlf_Der_p1_Inactive$FDR<0.2,]

# Save the DE results as CSV file
write.csv(all_results_qlf_Der_p1_Inactive, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Der_p1_Inactive.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_Der_p1_Inactive, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_qlf_Der_p1_Inactive.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_Der_p1_Inactive, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Inactive.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_Der_p1_Inactive,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_Der_p1_Inactive.rds"
)


###==== groupPoly_IC ====###

# Perform quasi-likelihood F test on each miRNA
qlf_groupPoly_IC <- glmQLFTest(fit, coef="groupPoly_IC")

# Calculate differential expression
all_results_qlf_groupPoly_IC <- as.data.frame(topTags(qlf_groupPoly_IC, n=nrow(qlf_groupPoly_IC)))

# Extract significant results 
sig_results_FDR005_qlf_groupPoly_IC <- all_results_qlf_groupPoly_IC[all_results_qlf_groupPoly_IC$FDR<0.05,]
sig_results_FDR02_qlf_groupPoly_IC <- all_results_qlf_groupPoly_IC[all_results_qlf_groupPoly_IC$FDR<0.2,]

# Save the DE results as CSV file
write.csv(all_results_qlf_groupPoly_IC, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_groupPoly_IC.csv", row.names = FALSE)
write.csv(sig_results_FDR005_qlf_groupPoly_IC, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR005_qlf_groupPoly_IC.csv", row.names = FALSE)
write.csv(sig_results_FDR02_qlf_groupPoly_IC, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/sig_results_FDR02_qlf_groupPoly_IC.csv", row.names = FALSE)

# Save DE object for MA plot
saveRDS(
  all_results_qlf_groupPoly_IC,
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/02_HDM_Dog/results/analysis_pipeline/run_01/07_DE_analysis/all_results_qlf_groupPoly_IC.rds"
)


