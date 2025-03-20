# Load necessary libraries
library(tidyverse)
library(RColorBrewer)
library(readxl)

# Read the GO enrichment data for each cluster
GO_cluster_ES <- read.table("Desktop/GO_ES_20250320.txt",header = T, fill = TRUE, sep = "\t")
GO_cluster_MS <- read.table("Desktop/GO_MS_20250320.txt",header = T, fill = TRUE, sep = "\t")

# Define pathways for each cluster
Pathway_cluster_ES <- c("macromolecule modification",
                      "post-translational protein modification",
                      "protein modification process",
                      "biological regulation")

Pathway_cluster_MS <- c("regulation of biological process",
                      "regulation of cellular process",
                      "biological regulation",
                      "regulation of transcription")

# Define cluster names and corresponding colors
cluster_names <- c(
  "ES",
  "MS"
)

celltype_colors <- c(
  "#8CA3C2"  ,
  "#BD746C"
)

# Function to create and save barplot for each cluster
create_barplot <- function(GO_data, pathways, cluster_name, fixed_max_limit) {
  # Filter for pathways of interest
  GO_data_filtered <- GO_data %>%
    filter(Term %in% pathways) %>%
    mutate(log_p_value = -log10(pvalue))  # Calculate -log10(p_value) and add as a new column
  
  # Print the updated data with the new column (log_p_value)
  print(head(GO_data_filtered))  # You can remove this line after checking
  
  # Calculate a fixed left margin (distance to the left edge of the bars)
  left_margin <- 0.1  # You can adjust this margin value
  
  # Get the color corresponding to the cluster_name
  cluster_color <- celltype_colors[which(cluster_names == cluster_name)]
  
  # Create the plot
  p <- ggplot(GO_data_filtered, aes(x = log_p_value, y = reorder(Term, log_p_value))) +
    geom_bar(stat = "identity", fill = cluster_color, width = 0.9) +  # Set the color dynamically
    # Apply dynamic x-axis for text, adjusting with a fixed left margin
    geom_text(aes(label = Term, x = left_margin), 
              hjust = 0, color = "black", size = 12) +  # Left-align text, dynamically set x position
    theme_classic() +
    theme(axis.title = element_blank(),   # Remove axis titles
          axis.text = element_blank(),    # Remove axis text
          axis.ticks = element_blank(),   # Remove axis ticks
          plot.title = element_blank(),   # Remove plot title
          panel.grid = element_blank(),   # Remove grid lines
          panel.background = element_blank(), # Remove background
          axis.line = element_blank()) +    # Ensure axis lines are removed
    scale_x_continuous(
      breaks = seq(0, fixed_max_limit, by = 5),    # Set x-axis ticks by 5 units
      limits = c(0, fixed_max_limit),              # Set fixed x-axis limits
      labels = scales::number_format()       # Display numeric values on x-axis
    ) +
    # Adding x-axis labels and the title for better readability
    labs(x = "-log10(p-value)", y = "Pathway") +
    theme(axis.text.x = element_text(size = 12, color = "black"))  # Display x-axis coordinates
  
  # Save the plot to a PDF file
  output_file <- paste0("Documents/GitHub/repli-ATAC-seq/Figures/", cluster_name, ".pdf")
  ggsave(output_file, plot = p, width = 20, height = 6)  # Adjust size as necessary
}

# Calculate the maximum log_p_value across all clusters
all_GO_data <- bind_rows(
  GO_cluster_ES %>% mutate(cluster = "ES"),
  GO_cluster_MS %>% mutate(cluster = "MS"),
)

# Calculate the maximum -log10(p-value) from all clusters and selected pathways
all_GO_data_filtered <- all_GO_data %>%
  filter(Term %in% c(Pathway_cluster_ES, Pathway_cluster_MS)) %>%
  mutate(log_p_value = -log10(pvalue))

fixed_max_limit <- max(all_GO_data_filtered$log_p_value, na.rm = TRUE) + 5  # Add some margin to the max value

# Cluster 2 barplot
create_barplot(GO_cluster_ES, Pathway_cluster_ES, "ES", fixed_max_limit)

# Cluster 3 barplot
create_barplot(GO_cluster_MS, Pathway_cluster_MS, "MS", fixed_max_limit)
