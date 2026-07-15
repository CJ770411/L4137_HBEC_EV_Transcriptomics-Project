# Load necessary libraries
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(plotly)
library(htmlwidgets)

#==========================#
# User Configuration 
#==========================#

logFC_threshold <- 1 # miRNA expression threshold for ORA

#==========================#
# Load DE object     
#==========================#


all_results_dose24 <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/all_results_dose24.rds"
)

#==========================#
# Volcano Plot 
#==========================#

all_results_dose24$significance <- "Not significant"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC > logFC_threshold] <- "Upregulated"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC < -logFC_threshold] <- "Downregulated"
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
  geom_vline(xintercept = c(1, -1), col = 'black', linetype = 'dashed', alpha = 0.6) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white')) +
  labs(x = "Log Fold Change", y = "-Log10 FDR") +
  scale_color_manual(values = c('Not signficant' = 'lightgrey', 'Upregulated' = '#E64B35', 'Downregulated' = '#0072B2'),
                     name = 'Significance') +
  coord_cartesian(xlim = c(x_max, x_min))

# Visualise volcano plot
print(volcano_plot)

# Save volcano plot
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/panels/volcano_plot.svg", plot = volcano_plot, width=10, height=7, units='in')

# Visualise interactive volcano plot
volcano_plot %>%
  ggplotly(tooltip = 'text')

# Save interactive volcano plot
saveWidget(
  ggplotly(volcano_plot, tooltip = "text"),
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_3/panels/volcano_plot.html",
  selfcontained = TRUE
)


