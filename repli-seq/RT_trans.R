# Set working directory if needed
# setwd("/your/working/dir")
setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation")

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

# Step 1: Load classified CSV
df <- read.csv("repli_phase_classified.csv")

# Step 2: 准备 WT 和各突变体的 RT 分类列
df_rt <- df %>%
  select(chr, start, end,
         wt = WT_phase,
         cr_sol1_5 = sol1_5_phase,
         cr_sol1_8 = sol1_8_phase,
         cr_tcx2_1 = tcx2_1_phase,
         cr_tcx2_3 = tcx2_3_phase) %>%
  mutate(across(everything(), ~ gsub("S", "", .))) %>%
  drop_na()

# Step 3: Stacked bar plot of WT-to-mutant transitions
df_long <- df_rt %>%
  pivot_longer(cols = starts_with("cr_"), names_to = "mutant", values_to = "mutant_RT")

transition_counts <- df_long %>%
  group_by(mutant, wt, mutant_RT) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(mutant, wt) %>%
  mutate(proportion = count / sum(count))

ggplot(transition_counts, aes(x = wt, y = proportion, fill = mutant_RT)) +
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Step 4: Define RT order and colors
rt_order <- c("E", "EM", "M", "ML", "L", "EL", "EML","Non-replication")
rt_colors <- c(
  "E" = "#2C5F9E", "EM" = "#68A0D8", "M" = "#95BE6C",
  "ML" = "#E4B660", "L" = "#E68364", "EL" = "#B784A7", "EML" = "#9B7EB3",
  "Non-replication" = "gray"
)

# Step 5: Alluvial plot function
plot_alluvial_pair <- function(mutant_col, mutant_label) {
  df_sub <- df_rt %>%
    select(wt, mutant = all_of(mutant_col)) %>%
    mutate(
      change_status = if_else(wt == mutant, "unchanged", "changed"),
      wt = factor(wt, levels = rt_order),
      mutant = factor(mutant, levels = rt_order)
    ) %>%
    filter(change_status == "changed") %>%
    group_by(wt, mutant) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(id = row_number())
  
  df_long <- df_sub %>%
    pivot_longer(cols = c("wt", "mutant"), names_to = "condition", values_to = "RT") %>%
    arrange(id, condition) %>%
    mutate(
      RT = factor(RT, levels = rt_order),
      alluvium = rep(df_sub$id, each = 2),
      Freq = rep(df_sub$Freq, each = 2),
      RT_wt = rep(df_sub$wt, each = 2)
    )
  
  ggplot(df_long,
         aes(x = condition, stratum = RT, alluvium = alluvium, y = Freq,
             fill = RT, label = RT)) +
    geom_flow(
      aes(fill = RT_wt),
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

# Step 6: Plot for each mutant
plot_alluvial_pair("cr_sol1_5", "cr_sol1_5")
plot_alluvial_pair("cr_sol1_8", "cr_sol1_8")
plot_alluvial_pair("cr_tcx2_1", "cr_tcx2_1")
plot_alluvial_pair("cr_tcx2_3", "cr_tcx2_3")

