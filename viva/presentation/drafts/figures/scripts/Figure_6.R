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

# Biologically plausible pathways
pathways_of_interest_up <- c(
  "Interleukin-4 and Interleukin-13 signaling",
  "Interleukin-17 signaling",
  "Signaling by Interleukins",
  "Signaling by TGFB family members",
  "Signaling by TGF-beta Receptor Complex",
  "TGF-beta receptor signaling activates SMADs",
  "Signaling by WNT",
  "Signaling by NOTCH",
  "Toll Like Receptor TLR1:TLR2 Cascade",
  "Toll Like Receptor 2 (TLR2) Cascade",
  "Toll Like Receptor TLR6:TLR2 Cascade",
  "Toll Like Receptor 5 (TLR5) Cascade",
  "Toll Like Receptor 10 (TLR10) Cascade",
  "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation"
)

# Duplicate pathway results to later create subset of biologically plausible pathways
epath_up_subset <- epath_up

# Create subset of biologically plausible pathways
epath_up_subset@result <- epath_up@result %>%
  filter(Description %in% pathways_of_interest_up)

# Create dotplot
dotplot_up <- dotplot(epath_up_subset, showCategory = 15) +
  labs_pubr() +
  theme_cleveland() +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = 'P.adj',
       size = 'Gene Count',
       title = 'Upregulated') +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 28),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white'),
        axis.text.y = element_text(size = 22)) 


# Visualise dotplot
print(dotplot_up)


#==========================#
# Panel B - Dotplot (downregulated)
#==========================#

epath_down <- readRDS(
  "~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/08_FE_indirect/epath_down_results.rds")

# Biologically plausible pathways
pathways_of_interest_down <- c(
  # Interleukin signalling
  "Interleukin-4 and Interleukin-13 signaling",
  "Interleukin-6 signaling",
  "Interleukin-7 signaling",
  "Interleukin-12 family signaling",
  "Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation",
  "Interleukin-17 signaling",
  "Signaling by Interleukins",
  
  # TGF-beta signalling
  "Signaling by TGFB family members",
  "Signaling by TGF-beta Receptor Complex",
  "TGF-beta receptor signaling activates SMADs",
  
  # Toll-like receptor signalling
  "Toll-like Receptor Cascades",
  "TRIF (TICAM1)-mediated TLR4 signaling",
  "Toll Like Receptor TLR1:TLR2 Cascade",
  "Toll Like Receptor 2 (TLR2) Cascade",
  "Toll Like Receptor TLR6:TLR2 Cascade",
  "Toll Like Receptor 3 (TLR3) Cascade",
  "Toll Like Receptor 4 (TLR4) Cascade",
  "Toll Like Receptor 5 (TLR5) Cascade",
  "Toll Like Receptor 7/8 (TLR7/8) Cascade",
  "Toll Like Receptor 9 (TLR9) Cascade",
  "Toll Like Receptor 10 (TLR10) Cascade",
  "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation",
  
  # Fc receptor signalling
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  "FCGR3A-mediated phagocytosis",
  "Fc epsilon receptor (FCERI) signaling",
  "FCERI mediated MAPK activation",
  
  # Pyroptosis
  "Defective pyroptosis",
  
  # Rho GTPase signalling
  "RHO GTPase cycle",
  "RHO GTPase Effectors",
  "RHO GTPases Activate Formins",
  "RHO GTPases Activate WASPs and WAVEs",
  "RHO GTPases activate PAKs",
  "RHO GTPases activate PKNs"
)

# Duplicate pathway results to later create subset of biologically plausible pathways
epath_down_subset <- epath_down

# Create subset of biologically plausible pathways
epath_down_subset@result <- epath_down@result %>%
  filter(Description %in% pathways_of_interest_down)


# Create dotplot
dotplot_down <- dotplot(epath_down_subset, showCategory = 15) +
  labs_pubr() +
  theme_cleveland() +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = 'P.adj',
       size = 'Gene Count',
       title = 'Downregulated') +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  theme(axis.ticks.length = unit(0.3, 'cm'),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 28),
        axis.title = element_text(size=20),
        legend.key.size = unit(1.0, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(fill = 'white'),
        axis.text.y = element_text(size = 22)) 


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
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/viva/presentation/drafts/figures/figure_6.pdf", plot = figure_6, height=8, width=12, dpi=600, units="in", scale=2)

