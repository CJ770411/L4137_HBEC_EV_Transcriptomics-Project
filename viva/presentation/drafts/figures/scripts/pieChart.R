library(RColorBrewer)
library(viridis)

# Create pie chart colour palette
pie_colours <- brewer.pal(6, "Dark2")


#==========================#
# Panel A - Pie chart (upregulated)
#==========================#

# Define each category (PubMed search term)
terms_125 <- c("Airway (5)", "Immune (85)", "Allergic (3)", "Inflammation (98)", "Epithelial (44)", "Extracellular Vesicle (89)")

# Define number of publications in each category
publications_125 <- c(5, 85, 3, 98, 44, 89)

# Prepare plot for saving
pdf("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/viva/presentation/drafts/figures/piechart_125.pdf", width = 7, height = 7)

# Create pie chart
pie(publications_125, labels = terms_125, col = pie_colours, density = 75, clockwise = TRUE)

# Save plot
dev.off()

#==========================#
# Panel B - Pie chart (upregulated)
#==========================#

# Define each category (PubMed search term)
terms_574 <- c("Airway (0)", "Immune (12)", "Allergic (1)", "Inflammation (19)", "Epithelial (13)", "Extracellular Vesicle (20)")

# Define number of publications in each category
publications_574 <- c(0, 12, 1, 19, 13, 20)

# Prepare plot for saving
pdf("~/GIT/L4137_HBEC_EV_Transcriptomics-Project/viva/presentation/drafts/figures/piechart_574.pdf", width = 7, height = 7)

# Create pie chart
pie(publications_574, labels = terms_574, col = pie_colours, density = 75, clockwise = TRUE)

# Save plot
dev.off()




