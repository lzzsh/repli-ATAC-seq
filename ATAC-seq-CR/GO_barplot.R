# ===============================
# 1️⃣ 加载依赖
# ===============================
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readxl)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

# ===============================
# 2️⃣ 读取你的 GO 富集文件
go_data <- read.table("Documents/GitHub/repli-ATAC-seq/data/GO_enrichment.txt",
                      header = TRUE, sep = "\t", fill = TRUE, quote = "")

# ===============================
# 3️⃣ 数据清洗与预处理
# ===============================
go_data <- go_data %>%
  rename(
    GO_ID = GO.term,
    Description = Description,
    P_value = P.value,
    FDR = FDR,
    Group = Group,
    Genes = GeneID
  ) %>%
  mutate(
    logP = -log10(as.numeric(P_value)),
    Group = case_when(
      Group == "biological_process" ~ "Biological Process",
      Group == "molecular_function" ~ "Molecular Function",
      Group == "cellular_component" ~ "Cellular Component",
      TRUE ~ Group
    )
  )

# ===============================
# 4️⃣ 过滤显著 GO 条目 (FDR < 0.05)
# ===============================
go_sig <- go_data %>%
  filter(FDR < 0.05)

# ===============================
# 5️⃣ 按每个 Group 分组，取前若干条
# ===============================
top_n_terms <- 10  # 每个组取前 10 个
go_top <- go_sig %>%
  group_by(Group) %>%
  slice_max(order_by = logP, n = top_n_terms)

# ===============================
# 6️⃣ 设定分组颜色
# ===============================
group_colors <- c(
  "Biological Process" = "#8CA3C2",
  "Molecular Function" = "#A58C7F",
  "Cellular Component" = "#BD746C"
)

# ===============================
# 7️⃣ 绘制 GO 富集条形图
# ===============================
p <- ggplot(go_top, aes(x = logP, y = reorder(Description, logP), fill = Group)) +
  geom_col(width = 0.8) +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
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

# ===============================
# 8️⃣ 保存结果
# ===============================
ggsave("Documents/GitHub/repli-ATAC-seq/Figures/GO_enrichment_barplot.pdf",
       plot = p, width = 8, height = 10, useDingbats = FALSE)

print("✅ 绘图完成：已输出 GO 富集条形图 PDF 文件！")