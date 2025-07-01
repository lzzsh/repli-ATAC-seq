# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation")

# Step 1: Read the classified replication phase results
df <- read.csv("repli_phase_classified.csv")

# Step 2: Extract the five sample columns and compute region length
df_rt <- df %>%
  dplyr::select(chr, start, end,
                WT = WT_phase,
                sol1_5 = sol1_5_phase,
                sol1_8 = sol1_8_phase,
                tcx2_1 = tcx2_1_phase,
                tcx2_3 = tcx2_3_phase) %>%
  mutate(length = end - start + 1)

# Step 3: Convert to long format and remove Non-replication and unknown
df_long <- df_rt %>%
  pivot_longer(cols = WT:tcx2_3, names_to = "sample", values_to = "RT") %>%
  filter(!RT %in% c("Non-replication", "unknown")) %>%
  mutate(
    RT = gsub("S", "", RT),  # Remove "S" from phase labels (e.g., "ES" → "E")
    RT = factor(RT, levels = c("E", "EM", "M", "ML", "L", "EL", "EML"))  # Set factor levels for order
  )

# Step 3.5: Summarize number and total length of Non-replication regions per sample
nonrep_summary <- df_rt %>%
  pivot_longer(cols = WT:tcx2_3, names_to = "sample", values_to = "RT") %>%
  filter(RT == "Non-replication") %>%
  group_by(sample) %>%
  summarise(
    region_count = n(),              # number of Non-replication regions
    total_length = sum(length),     # total genomic length of Non-replication regions
    .groups = "drop"
  )

# Print summary of Non-replication per sample
print(nonrep_summary)

# Step 4: Calculate total RT segment length and proportion per sample (excluding Non-replication)
pie_data <- df_long %>%
  group_by(sample, RT) %>%
  summarise(total_length = sum(length), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(
    proportion = total_length / sum(total_length),
    RT = fct_rev(RT)  # Reverse order for consistent pie chart layering
  )

# Step 5: Define custom color palette for RT phases
rt_colors <- c(
  "E" = "#2C5F9E", "EM" = "#68A0D8", "M" = "#95BE6C",
  "ML" = "#E4B660", "L" = "#E68364", "EL" = "#B784A7", "EML" = "#9B7EB3"
)

# Step 6: Draw multiple pie charts using facet_wrap (one per sample)
pie_plot <- ggplot(pie_data, aes(x = "", y = proportion, fill = RT)) +
  geom_bar(width = 1, stat = "identity", color = "black") +  # pie sector with borders
  coord_polar(theta = "y") +                                 # convert bar to pie chart
  facet_wrap(~ sample, ncol = 3) +                           # one pie per sample
  scale_fill_manual(values = rt_colors) +                    # apply color scheme
  theme_void(base_size = 14) +                               # clean background
  theme(
    strip.text = element_text(face = "bold", size = 14),     # facet titles
    legend.position = "right"                                # show legend on right
  ) +
  labs(
    title = "Replication Timing Segment Composition",
    fill = "RT Phase"
  )

# Display pie charts
print(pie_plot)

# Optional: Save pie plot to PDF
# ggsave("RT_phase_pie_by_sample.pdf", pie_plot, width = 10, height = 6)