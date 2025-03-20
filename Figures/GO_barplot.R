# Load necessary libraries
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(pheatmap)

# Read enrichment data
enrich <- read_excel("Documents/GitHub/repli-ATAC-seq/data/enrich.xlsx")

# Determine the main stage based on the highest odds ratio, including ESMS
enrich <- enrich %>%
  mutate(
    main_stage = case_when(
      odds_ratio_ES >= odds_ratio_ESMS & odds_ratio_ES >= odds_ratio_MS & odds_ratio_ES >= odds_ratio_LS ~ "ES",
      odds_ratio_ESMS >= odds_ratio_ES & odds_ratio_ESMS >= odds_ratio_MS & odds_ratio_ESMS >= odds_ratio_LS ~ "ESMS",
      odds_ratio_MS >= odds_ratio_ES & odds_ratio_MS >= odds_ratio_ESMS & odds_ratio_MS >= odds_ratio_LS ~ "MS",
      odds_ratio_LS >= odds_ratio_ES & odds_ratio_LS >= odds_ratio_ESMS & odds_ratio_LS >= odds_ratio_MS ~ "LS"
    )
  )

# Sort the data by stage (ES > ESMS > MS > LS) and descending odds ratio
enrich <- enrich %>%
  arrange(
    factor(main_stage, levels = c("ES", "ESMS", "MS", "LS")),  # Ensure ESMS is correctly positioned
    desc(if_else(main_stage == "ES", odds_ratio_ES,
                 if_else(main_stage == "ESMS", odds_ratio_ESMS,
                         if_else(main_stage == "MS", odds_ratio_MS, odds_ratio_LS))))
  )

# Extract relevant columns for heatmap
heatmap_data <- enrich %>%
  select(odds_ratio_ES, odds_ratio_ESMS, odds_ratio_MS, odds_ratio_LS) %>%
  as.matrix()

# Set row names as transcription factor (TF) names
rownames(heatmap_data) <- enrich$TF

# Plot heatmap
p <- pheatmap(
  heatmap_data,
  cluster_rows = FALSE,  # Maintain the sorted order
  cluster_cols = FALSE,  # Keep original column order
  show_rownames = FALSE, 
  cellwidth = 100,  # Adjust cell width
  cellheight = 0.008,
  scale = "row",  # Normalize row-wise
  main = "Stage-Specific Enrichment of eCAAS"
)

# Save heatmap as PDF
ggsave("Documents/GitHub/repli-ATAC-seq/Figures/enrichment_eCAAS.pdf", p, width = 8, height = 8, useDingbats = FALSE)

# Read GO enrichment data for each cluster
GO_cluster_ES <- read.table("Desktop/GO_ES_20250320.txt", header = TRUE, fill = TRUE, sep = "\t")
GO_cluster_ESMS <- read.table("Desktop/GO_ESMS_20250320.txt", header = TRUE, fill = TRUE, sep = "\t")
GO_cluster_MS <- read.table("Desktop/GO_MS_20250320.txt", header = TRUE, fill = TRUE, sep = "\t")

# Define pathways for each cluster
Pathway_cluster_ES <- c("macromolecule modification",
                        "post-translational protein modification",
                        "protein modification process",
                        "biological regulation")

Pathway_cluster_ESMS <- c("cell communication",
                          "regulation of signal transduction",
                          "biological regulation",
                          "regulation of transcription")

Pathway_cluster_MS <- c("regulation of biological process",
                        "regulation of cellular process",
                        "biological regulation",
                        "regulation of transcription")

# Define cluster names and corresponding colors
cluster_names <- c("ES", "ESMS", "MS")

celltype_colors <- c("#8CA3C2", "#A58C7F", "#BD746C")  # Assign colors for ES, ESMS, and MS

# Function to create and save a bar plot for each cluster
create_barplot <- function(GO_data, pathways, cluster_name, fixed_max_limit) {
  # Filter for pathways of interest
  GO_data_filtered <- GO_data %>%
    filter(Term %in% pathways) %>%
    mutate(log_p_value = -log10(pvalue))  # Convert p-value to -log10 scale
  
  # Print filtered data for verification (can be removed)
  print(head(GO_data_filtered))
  
  # Define a fixed left margin for text alignment
  left_margin <- 0.1
  
  # Get the corresponding color for the cluster
  cluster_color <- celltype_colors[which(cluster_names == cluster_name)]
  
  # Create the bar plot
  p <- ggplot(GO_data_filtered, aes(x = log_p_value, y = reorder(Term, log_p_value))) +
    geom_bar(stat = "identity", fill = cluster_color, width = 0.9) +  # Set bar color
    geom_text(aes(label = Term, x = left_margin), 
              hjust = 0, color = "black", size = 12) +  # Left-align text dynamically
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()) +
    scale_x_continuous(
      breaks = seq(0, fixed_max_limit, by = 5),  # Set x-axis ticks
      limits = c(0, fixed_max_limit),            # Fix x-axis limits
      labels = scales::number_format()           # Format numbers
    ) +
    labs(x = "-log10(p-value)", y = "Pathway") +
    theme(axis.text.x = element_text(size = 12, color = "black"))
  
  # Save the plot as a PDF
  output_file <- paste0("Documents/GitHub/repli-ATAC-seq/Figures/", cluster_name, ".pdf")
  ggsave(output_file, plot = p, width = 20, height = 6)
}

# Combine all GO data for global scaling of x-axis
all_GO_data <- bind_rows(
  GO_cluster_ES %>% mutate(cluster = "ES"),
  GO_cluster_ESMS %>% mutate(cluster = "ESMS"),
  GO_cluster_MS %>% mutate(cluster = "MS")
)

# Calculate the maximum -log10(p-value) for all clusters
all_GO_data_filtered <- all_GO_data %>%
  filter(Term %in% c(Pathway_cluster_ES, Pathway_cluster_ESMS, Pathway_cluster_MS)) %>%
  mutate(log_p_value = -log10(pvalue))

fixed_max_limit <- max(all_GO_data_filtered$log_p_value, na.rm = TRUE) + 5  # Add margin to max value

# Generate bar plots for all clusters
create_barplot(GO_cluster_ES, Pathway_cluster_ES, "ES", fixed_max_limit)
create_barplot(GO_cluster_ESMS, Pathway_cluster_ESMS, "ESMS", fixed_max_limit)
create_barplot(GO_cluster_MS, Pathway_cluster_MS, "MS", fixed_max_limit)