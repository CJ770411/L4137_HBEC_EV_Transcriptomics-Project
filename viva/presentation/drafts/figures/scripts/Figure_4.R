# Load necessary libraries
library(readr)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(patchwork)

#==========================#
# Panel C - Differential Expression Barchart
#==========================#


# Generate Bar Chart Data (Protease)   

# Load primary HBEC DE data
dose0_005_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose0005.csv")
dose0_08_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose008.csv")
dose2_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose2.csv")
dose24_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/01_Protease/07_DE_analysis/sig_results_FDR005_dose24.csv")

# Add condition labels
dose0_005_data <- dose0_005_data %>% 
  count() %>%
  mutate(Condition = 'Dose 0.005')

dose0_08_data <- dose0_08_data %>% 
  count() %>%
  mutate(Condition = 'Dose 0.08')

dose2_data <- dose2_data %>% 
  count() %>%
  mutate(Condition = 'Dose 2')

dose24_data <- dose24_data %>% 
  count() %>%
  mutate(Condition = 'Dose 24')



# Generate Bar Chart Data (HDM_Dog)   

# Load DE data
can_f_1_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Can_F_1.csv")
der_p1_active_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Active.csv")
der_p1_inactive_data <- read_csv("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/write_up/drafts/results/02_HDM_Dog/07_DE_analysis/sig_results_FDR02_qlf_Der_p1_Inactive.csv")

# Add condition labels
can_f_1_data <- can_f_1_data %>% 
  count() %>%
  mutate(Condition = 'Can f 1')

der_p1_active_data <- der_p1_active_data %>% 
  count() %>%
  mutate(Condition = 'Der P1 Active')

der_p1_inactive_data <- der_p1_inactive_data %>% 
  count() %>%
  mutate(Condition = 'Der P1 Inactive')



# Create Bar Chart 

# Combine condition data
boxplot_data <- rbind(dose0_005_data,
                      dose0_08_data,
                      dose2_data,
                      dose24_data,
                      can_f_1_data,
                      der_p1_active_data,
                      der_p1_inactive_data)

# Preserve condition order
boxplot_data$Condition <- factor(
  boxplot_data$Condition,
  levels = c(
    "Dose 0.005",
    "Dose 0.08",
    "Dose 2",
    "Dose 24",
    "Can f 1",
    "Der P1 Active",
    "Der P1 Inactive"
  )
)

# Set column name as DE miRNA count
colnames(boxplot_data)[colnames(boxplot_data) == 'n'] <- 'Count'

# Plot bar chart
diffexpr_barchart <- ggplot(boxplot_data, aes(x = Condition, y = Count, fill = Condition)) +
  geom_col() +
  labs(x = 'Condition',
       y = 'DE miRNA Count') +
  theme_classic() + 
  labs_pubr() +
  scale_fill_manual(values = c(
    'Dose 0.005'      = '#D1EEE8',
    'Dose 0.08'       = '#7BC8B6',
    'Dose 2'          = '#2A9D8F',
    'Dose 24'         = '#146C63',
    'Can f 1'         = '#E8D5F2',
    'Der P1 Active'   = '#C7A0D8',
    'Der P1 Inactive' = '#9C6BB3',
    'Poly I:C'        = '#6A3D9A'
  )) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 38),
        axis.title.y = element_text(size = 42),
        axis.title.x = element_text(size = 42)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_y_continuous(
    breaks = seq(0, 30, by = 5))

# Visualise bar chart
print(diffexpr_barchart)


# Save figure
ggsave("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/viva/presentation/drafts/figures/figure_4.pdf", plot = diffexpr_barchart, height=6, width=16, dpi=600, units="in", scale=2)


