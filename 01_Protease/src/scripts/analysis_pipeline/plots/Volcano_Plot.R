# Load necessary libraries
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(plotly)
library(htmlwidgets)



#==========================#
# Load DE object     
#==========================#


all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/07_DE_analysis/all_results_dose24.rds"
)

#==========================#
# Volcano Plot 
#==========================#

all_results_dose24$significance <- "Not significant"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC > 0] <- "Upregulated"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC < 0] <- "Downregulated"
all_results_dose24$significance <- factor(all_results_dose24$significance, levels = c("Upregulated", "Downregulated", "Not significant"))



# Set axis limits
x_max <- ceiling(max(abs(all_results_dose24$logFC)))
x_min <- -x_max

volcano_plot <- ggplot(all_results_dose24, aes(x = logFC, y = -log10(FDR), color = significance, text = paste('Symbol:', Symbol,
                                                                                                         '<br>logFC:', round(logFC, 3),
                                                                                                         '<br>FDR:', signif(FDR, 3)))) +
  geom_point(size = 4) +
  theme_classic() +
  labs_pubr() +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed', alpha = 0.6) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white')) +
  labs(x = "Log Fold Change", y = "-Log10 FDR") +
  scale_color_manual(values = c('Not signficant' = 'lightgrey', 'Upregulated' = '#E64B35', 'Downregulated' = '#4DBBD5'),
                     name = 'Significance') +
  coord_cartesian(xlim = c(x_max, x_min))

# Visualise volcano plot
print(volcano_plot)

# Save volcano plot
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/plots/volcano_plot.svg", plot = volcano_plot, width=10, height=7, units='in')

# Visualise interactive volcano plot
volcano_plot %>%
  ggplotly(tooltip = 'text')

# Save interactive volcano plot
saveWidget(
  ggplotly(volcano_plot, tooltip = "text"),
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/plots/volcano_plot.html",
  selfcontained = TRUE
)


