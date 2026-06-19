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
# MA Plot 
#==========================#

all_results_dose24$significance <- "Not significant"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC > 0] <- "Upregulated"
all_results_dose24$significance[all_results_dose24$FDR < 0.05 & all_results_dose24$logFC < 0] <- "Downregulated"
all_results_dose24$significance <- factor(all_results_dose24$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

# Set axis limits
x_min <- floor(min(all_results_dose24$logCPM))
x_max <- ceiling(max(all_results_dose24$logCPM))

# Create MA plot
ma_plot <- ggplot(all_results_dose24, aes(x = logCPM, y = logFC, color = significance, size = significance, alpha = significance,
                               text = paste(
                                 "Symbol:", Symbol,
                                 "<br>logFC:", round(logFC, 3),
                                 "<br>logCPM:", round(logCPM, 3),
                                 "<br>FDR:", signif(FDR, 3)
                                 ))) +
  geom_point() +
  scale_color_manual(
    values = c("Not significant" = "gray80",
               "Downregulated" = "#0072B2",
               "Upregulated" = "#D55E00")
  ) +
  scale_size_manual(    
    values = c("Not significant" = 3,
               "Downregulated" = 4,
               "Upregulated" = 4)) +
  scale_alpha_manual(    
    values = c("Not significant" = 0.4,
               "Downregulated" = 0.8,
               "Upregulated" = 0.8)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  labs(x = "Average Log CPM",
       y = "Log Fold Change",
       color = '',
       size = '',
       alpha = '') +
  theme_pubclean() +
  labs_pubr() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_continuous(
    breaks = seq(x_min, x_max, by = 1))

# Visualise MA plot
print(ma_plot)


# Save MA plot
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/plots/MA_plot.svg", plot = ma_plot, width=10, height=7, units='in')

# Visualise interactive MA plot
ma_plot %>%
  ggplotly(tooltip = "text")

# Save interactive MA plot
saveWidget(
  ggplotly(ma_plot, tooltip = "text"),
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/01_Protease/results/analysis_pipeline/run_01/plots/MA_plot.html",
  selfcontained = TRUE
)

