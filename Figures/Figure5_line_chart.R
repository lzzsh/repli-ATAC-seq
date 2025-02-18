# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)

# Create the dataset
df <- data.frame(
  Phase = c("ES","ESMS","MS","MSLS","LS"),  # Replication phases
  wt = c(16307, 4434, 19633, 1101, 29596),  # Peak counts for wild-type (wt)
  xw11_OE = c(14041, 5240, 20606, 1466, 32022),  # Peak counts for xw11.OE
  xw11_CR = c(15239, 2055, 24943, 1364, 33113),  # Peak counts for xw11.CR
  x39_CR = c(16649, 4377, 20985, 1343, 33931)  # Peak counts for x39.CR
)

# Transform the dataset from wide format to long format
df_long <- df %>%
  pivot_longer(cols = -Phase, names_to = "Batch", values_to = "Count")

# Define the order of the batches
df_long$Batch <- factor(df_long$Batch, levels = c("wt", "xw11_OE", "xw11_CR", "x39_CR"))

# Define custom colors for each phase
custom_colors <- c("ES" = "#2C5F9E",  # Dark blue (Early phase)
                   "MS" = "#95BE6C",  # Green (Mid phase)
                   "LS" = "#E68364",
                   "MSLS" = "#E4B660",
                   "ESMS" = "#68A0D8")  # Orange-red (Late phase)

# Generate the line plot
p <- ggplot(df_long, aes(x = Batch, y = Count, group = Phase, color = Phase)) +
  geom_line(size = 2.5, linetype = "solid") +  # Adjust line thickness
  geom_point(size = 4, alpha = 0.9) +  # Add data points with transparency
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_classic() +
  labs(x = "Sample",  # X-axis label
       y = "Number of peaks",  # Y-axis label
       color = "Phase") +  # Legend title
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate X-axis labels
    axis.title.x = element_blank(),  # Remove X-axis title
    axis.title.y = element_text(size = 12),  # Adjust Y-axis title size
    axis.text.y = element_text(size = 12),  # Adjust Y-axis text size
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center align plot title
    legend.title = element_text(face = "bold", size = 12),  # Bold legend title
    legend.key.size = unit(1.5, "cm"),  # Adjust legend key size
    legend.justification = c(0.4, 0.5),  # Adjust legend alignment
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "top"  # Place legend at the top
  )

# Save the plot as a PDF file
ggsave("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/line_chart.pdf",
       p, width = 8, height = 6)


# # Load required libraries
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(ggalluvial)

# # Create the dataset
# df <- data.frame(
#   Phase = c("ES", "MS", "LS"),  # Replication phases
#   wt = c(16004, 17104, 24708),  # Peak counts for wild-type (wt)
#   xw11_OE = c(11061, 18225, 25686),  # Peak counts for xw11.OE
#   xw11_CR = c(11648, 21272, 26926),  # Peak counts for xw11.CR
#   x39_CR = c(13138, 17649, 27016)  # Peak counts for x39.CR
# )

# # Transform the dataset from wide format to long format
# df_long <- df %>%
#   pivot_longer(cols = -Phase, names_to = "Batch", values_to = "Count")

# # Define the order of the batches
# df_long$Batch <- factor(df_long$Batch, levels = c("wt", "xw11_OE", "xw11_CR", "x39_CR"))

# # Define custom colors for each phase
# custom_colors <- c("ES" = "#2C5F9E",  # Dark blue (Early phase)
#                    "MS" = "#95BE6C",  # Green (Mid phase)
#                    "LS" = "#E68364")  # Orange-red (Late phase)

# # Generate the line plot
# p <- ggplot(df_long, aes(x = Batch, y = Count, group = Phase, color = Phase)) +
#   geom_line(size = 2.5, linetype = "solid") +  # Adjust line thickness
#   geom_point(size = 4, alpha = 0.9) +  # Add data points with transparency
#   scale_color_manual(values = custom_colors) +  # Apply custom colors
#   theme_classic() +
#   labs(x = "Sample",  # X-axis label
#        y = "Number of peaks",  # Y-axis label
#        color = "Phase") +  # Legend title
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate X-axis labels
#     axis.title.x = element_blank(),  # Remove X-axis title
#     axis.title.y = element_text(size = 12),  # Adjust Y-axis title size
#     axis.text.y = element_text(size = 12),  # Adjust Y-axis text size
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center align plot title
#     legend.title = element_text(face = "bold", size = 12),  # Bold legend title
#     legend.key.size = unit(1.5, "cm"),  # Adjust legend key size
#     legend.justification = c(0.4, 0.5),  # Adjust legend alignment
#     legend.text = element_text(size = 12),  # Adjust legend text size
#     legend.position = "top"  # Place legend at the top
#   )

# # Save the plot as a PDF file
# ggsave("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/line_chart.pdf",
#        p, width = 8, height = 6)