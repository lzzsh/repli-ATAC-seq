library(dplyr)

# List of GFF file paths
gff_files <- list(
  wt = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_wt.gff3",
  cr_sol1_5 = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_sol1_5.gff3",
  cr_sol1_8 = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_sol1_8.gff3",
  cr_tcx2_1 = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_tcx2_1.gff3",
  cr_tcx2_3 = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_tcx2_3.gff3"
)

# Define RT segmentation parser and splitter function
process_gff <- function(path) {
  gff <- read.table(path, comment.char = "#")[, c(1, 4, 5, 9)]
  colnames(gff) <- c("chr", "start", "end", "RT")
  gff$RT <- sapply(gff$RT, function(x) unlist(strsplit(x, ";"))[2])
  gff$RT <- sapply(gff$RT, function(x) unlist(strsplit(x, "="))[2])
  gff <- gff[!(gff$chr %in% c("chrUn", "chrSy")), ]
  gff$start <- as.integer(gff$start)
  gff$end <- as.integer(gff$end)
  gff <- gff[!is.na(gff$RT), ]
  gff$RT <- gsub("S", "", gff$RT)
  
  # Split long intervals
  split_gff_row <- function(row) {
    len <- row$end - row$start + 1
    if (len <= 1000) return(data.frame(
      chr = row$chr,
      start = row$start,
      end = row$end,
      RT = row$RT,
      stringsAsFactors = FALSE
    ))
    n_parts <- len %/% 1000
    starts <- seq(row$start, by = 1000, length.out = n_parts)
    ends <- starts + 999
    data.frame(
      chr = row$chr,
      start = starts,
      end = ends,
      RT = row$RT,
      stringsAsFactors = FALSE
    )
  }
  
  gff_split <- bind_rows(lapply(1:nrow(gff), function(i) split_gff_row(gff[i, ])))
  gff_split <- gff_split %>% arrange(chr, start)
  return(gff_split)
}

# Apply the function to all files
gff_results <- lapply(gff_files, process_gff)

merged_df <- Reduce(function(x, y) full_join(x, y, by = c("chr", "start", "end")),
                    lapply(names(gff_results), function(name) {
                      g <- gff_results[[name]]
                      colnames(g)[4] <- name  # Rename RT column
                      return(g)
                    }))

library(tidyverse)

# Step 1: Remove NA rows
df <- merged_df %>% drop_na()

# Step 2: Convert to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("cr_"), names_to = "mutant", values_to = "mutant_RT")

# Step 3: Count transitions
transition_counts <- df_long %>%
  group_by(mutant, wt, mutant_RT) %>%
  summarise(count = n(), .groups = "drop")

# Step 4: Calculate proportions within each WT group for each mutant
transition_props <- transition_counts %>%
  group_by(mutant, wt) %>%
  mutate(proportion = count / sum(count))

# Step 5: Plot stacked barplot by proportion
ggplot(transition_props, aes(x = wt, y = proportion, fill = mutant_RT)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ mutant, ncol = 2) +
  labs(
    title = "RT Transition Proportions from WT to Mutants",
    x = "RT State in WT",
    y = "Proportion",
    fill = "RT in Mutant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

library(tidyverse)
library(ggalluvial)

# Step 1: Remove NA rows
df <- merged_df %>% drop_na()

# Step 2: Select relevant columns
df_use <- df %>% dplyr::select(wt, cr_sol1_5, cr_sol1_8, cr_tcx2_1, cr_tcx2_3)

# Step 3: Define RT display order
rt_order <- c("E", "EM", "M", "ML", "L", "EL", "EML")

# Step 4: Define plotting function
plot_alluvial_pair <- function(mutant_col, mutant_label) {
  df_sub <- df_use %>%
    dplyr::select(wt, mutant = all_of(mutant_col)) %>%
    mutate(
      change_status = if_else(wt == mutant, "unchanged", "changed"),
      wt = factor(wt, levels = rt_order),
      mutant = factor(mutant, levels = rt_order)
    ) %>%
    filter(change_status == "changed") %>%  # Only changed
    group_by(wt, mutant) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(id = row_number())
  
  # Convert to long format for ggalluvial
  df_long <- df_sub %>%
    pivot_longer(cols = c("wt", "mutant"), names_to = "condition", values_to = "RT") %>%
    arrange(id, condition) %>%
    mutate(
      RT = factor(RT, levels = rt_order),
      alluvium = rep(df_sub$id, each = 2),
      Freq = rep(df_sub$Freq, each = 2),
      RT_wt = rep(df_sub$wt, each = 2)  # ← WT-based color
    )
  
  # Custom color palette for RT states
  rt_colors <- c(
    "E" = "#2C5F9E",
    "EM" = "#68A0D8",
    "M" = "#95BE6C",
    "ML" = "#E4B660",
    "L" = "#E68364",
    "EL" = "#B784A7",
    "EML" = "#9B7EB3"
  )
  
  # Generate alluvial plot
  ggplot(df_long,
         aes(x = condition, stratum = RT, alluvium = alluvium, y = Freq,
             fill = RT, label = RT)) +
    geom_flow(
      aes(fill = RT_wt),  # ← Flow colored by WT RT
      stat = "alluvium",
      lode.guidance = "forward",
      alpha = 0.5,
      curve_type = "sigmoid"
    ) +
    geom_stratum(aes(fill = RT), alpha = 0.9, color = "black") +
    scale_fill_manual(values = rt_colors, drop = FALSE) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("RT Transitions (Changed Only): WT →", mutant_label),
      x = "Condition", y = "Count",
      fill = "RT State"
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
}

# Step 5: Generate plots for each mutant
plot_alluvial_pair("cr_sol1_5", "cr_sol1_5")
plot_alluvial_pair("cr_sol1_8", "cr_sol1_8")
plot_alluvial_pair("cr_tcx2_1", "cr_tcx2_1")
plot_alluvial_pair("cr_tcx2_3", "cr_tcx2_3")