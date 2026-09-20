library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(readr)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

go_data <- read.csv("./GO_SOL1-8_down.csv",
                      header = TRUE, fill = TRUE)

# 去除列名前后空格
colnames(go_data) <- trimws(colnames(go_data))

# 自动匹配关键列名
go_data <- go_data %>%
  rename_with(~ "GO_ID", grep("^GO", names(go_data), ignore.case = TRUE, value = TRUE)) %>%
  rename_with(~ "P_value", grep("^P(\\.|_)?value", names(go_data), ignore.case = TRUE, value = TRUE)) %>%
  rename_with(~ "FDR", grep("^FDR", names(go_data), ignore.case = TRUE, value = TRUE)) %>%
  rename_with(~ "Group", grep("^Group", names(go_data), ignore.case = TRUE, value = TRUE)) %>%
  rename_with(~ "Description", grep("^Desc", names(go_data), ignore.case = TRUE, value = TRUE)) %>%
  rename_with(~ "Genes", grep("^Gene", names(go_data), ignore.case = TRUE, value = TRUE))

go_data <- go_data %>%
  mutate(
    logP = -log10(as.numeric(P_value)),
    Group = case_when(
      Group == "biological_process" ~ "Biological Process",
      Group == "molecular_function" ~ "Molecular Function",
      Group == "cellular_component" ~ "Cellular Component",
      TRUE ~ Group
    )
  )

go_sig <- go_data %>%
  filter(FDR < 0.05)

top_n_terms <- 3  # 每个组取前 10 个
go_top <- go_sig %>%
  group_by(Group) %>%
  slice_max(order_by = logP, n = top_n_terms)

# 设定分组颜色
group_colors <- c(
  "Biological Process" = "#8CA3C2",
  "Molecular Function" = "#A58C7F",
  "Cellular Component" = "#BD746C"
)

# 绘制 GO 富集条形图
p <- ggplot(go_top, aes(x = logP, y = reorder(Description, logP), fill = Group)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = expression(-log[10]("P-value")),
    y = "GO term",
    title = "Top Enriched GO Terms by Category"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )


# 保存结果
ggsave("./sol1-8_GO_enrichment_barplot_down.pdf",
       plot = p, width = 8, height = 5, useDingbats = FALSE)
