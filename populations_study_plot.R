# Install and load necessary packages
library(UpSetR)
library(tidyverse)
library(tidyverse)
library(gridExtra)
library(grid)
library(svglite)

# Load the data
data <- read.delim("/Users/yperez/work/pgt-pangenome/blast_canonical-count-tables/1-2mismatches_peptides_match_info/peps_merged_populations.tsv", header = TRUE, sep = "\t")

# Ensure your data has the appropriate column names
# Assuming columns are named 'peptide_sequence', 'gene', and 'population_category'

# Prepare the data for UpSet plot
data_wide <- data %>%
  distinct(genes, population_group) %>%
  mutate(presence = 1) %>%
  spread(population_group, presence, fill = 0)


# Prepare the data for peptide sequences
data_peptides <- data %>%
  distinct(peptide, population_group) %>%
  mutate(presence = 1) %>%
  spread(population_group, presence, fill = 0)

# Prepare the data for genes
data_genes <- data %>%
  distinct(genes, population_group) %>%
  mutate(presence = 1) %>%
  spread(population_group, presence, fill = 0)

# Define a vector of colors for the population groups
population_colors <- c("#48D1CC", "#9370DB", "#1E90FF", "#4682B4", "#7a96ac", "#DAA520")


# Create an UpSet plot for peptide sequences
p1 <- UpSetR::upset(data_peptides, sets = colnames(data_peptides)[-1], main.bar.color = "#1f77b4",sets.bar.color = population_colors,order.by = "freq", empty.intersections = "on")
p2 <- UpSetR::upset(data_genes, sets = colnames(data_genes)[-1], main.bar.color = "#1f77b4", sets.bar.color = population_colors, order.by = "freq",empty.intersections = "on")

# Convert the UpSet plots to ggplot objects
p1_plot <- cowplot::plot_grid(NULL, p1$BoxPlots, p1$Sizes, ncol = 1, align = 'v', rel_heights = c(0.2, 0.5, 0.3))
p2_plot <- cowplot::plot_grid(NULL, p1$Main_bar, p1$Matrix, ncol = 1, align = 'v', rel_heights = c(0.2, 0.5, 0.3))
p3_plot <- cowplot::plot_grid(NULL, p2$BoxPlots, p2$Sizes, ncol = 1, align = 'v', rel_heights = c(0.2, 0.5, 0.3))
p4_plot <- cowplot::plot_grid(NULL, p2$Main_bar, p2$Matrix, ncol = 1, align = 'v', rel_heights = c(0.2, 0.5, 0.3))


# Arrange the plots side by side
final_plot <- cowplot::plot_grid(p1_plot, p2_plot, p3_plot, p4_plot, ncol = 4, rel_widths = c(0.3, 1.7, 0.3, 1.7))

# Save the final plot as PNG
ggsave("populations_variants_genes.png", final_plot, width = 12, height = 6, units = "in", dpi = 300)

# Save the final plot as SVG
ggsave("populations_variants_genes.svg", final_plot, width = 12, height = 6, units = "in")


