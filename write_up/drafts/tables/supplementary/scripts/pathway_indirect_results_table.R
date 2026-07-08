library(dplyr)

# Load indirect pathway analysis results table
epath_up_results_table <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/08_FE_indirect/epath_up_results_table.rds"
)

epath_down_results_table <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/08_FE_indirect/epath_down_results_table.rds"
)

# Filter results (Padj < 0.05; q < 0.02)
epath_up_results_table_filtered <- epath_up_results_table %>% 
  filter(p.adjust < 0.02,
         qvalue < 0.02)

epath_down_results_table_filtered <- epath_down_results_table %>% 
  filter(p.adjust < 0.02,
         qvalue < 0.02)


# Save results tables as CSV
write.csv(epath_up_results_table_filtered, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/tables/supplementary/pathway_indirect_results_table_up.csv")
write.csv(epath_down_results_table_filtered, "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/tables/supplementary/pathway_indirect_results_table_down.csv")
