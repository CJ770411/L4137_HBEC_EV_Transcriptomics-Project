# Load necessary libraries
library(patchwork)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)

#==========================#
# Panel A - Dotplot (upregulated)
#==========================#

# Load target gene-level enrichment results
epath_up <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/08_FE_indirect/epath_up_results.rds")

# Create dotplot
dotplot_up <- dotplot(epath_up, showCategory = 30) +
  labs_pubr() +
  theme_cleveland() +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = 'P.adj',
       size = 'Gene Count') +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white'),
        axis.text.y = element_text(size = 13)) 


# Visualise dotplot
print(dotplot_up)



#==========================#
# Panel B - Dotplot (downregulated)
#==========================#

epath_down <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/08_FE_indirect/epath_down_results.rds")


# Create dotplot
dotplot_down <- dotplot(epath_down, showCategory = 30) +
  labs_pubr() +
  theme_cleveland() +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = 'P.adj',
       size = 'Gene Count') +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white'),
        axis.text.y = element_text(size = 13)) 


# Visualise dotplot
print(dotplot_down)


#==========================#
# Combine Plot to Create Figure
#==========================#


# Combine plots
figure_6 <- (dotplot_up + dotplot_down) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.key.size = unit(1.2, "cm"),
        plot.tag = element_text(size = 22, face = 'bold'))

# Visualise figure
print(figure_6)

# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/figures/main/figure_6/figure_6.pdf", plot = figure_6, height=8, width=12, dpi=600, units="in", scale=2)

