############################################################
## 0. 加载依赖包
############################################################
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggalluvial)


############################################################
### 1. 读取 GFF3 并提取基因的 TSS 信息
############################################################

gtf_data <- import('/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3') %>% 
  as.data.frame()

gene_data <- gtf_data %>%
  filter(!is.na(gtf_data[,11]), gtf_data[,7] == "gene") %>%
  dplyr::select(
    chrom  = 1,
    start  = 2,
    end    = 3,
    strand = 5,
    Geneid = 11
  )

# 用 start 生成 1bp TSS 用于匹配
gene_data$end <- gene_data$start
gene_data$start <- gene_data$end - 1


############################################################
### 2. 读取 TSS 对应的复制时期 (RT)
############################################################

TSS_RT <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/TSS_RT/TSS_RT_bin.bed") %>%
  dplyr::select(
    chrom   = V1,
    start   = V2,
    end     = V3,
    strand  = V6,
    RT_WT   = V7,
    RT_SOL1 = V9,
    RT_TCX2 = V11
  )


############################################################
### 3. 基因与 RT 进行匹配
############################################################

gene_RT_full <- TSS_RT %>%
  inner_join(gene_data, by = c("chrom","start","end","strand")) %>%
  distinct(Geneid, RT_WT, RT_TCX2)


############################################################
### 4. 引入 DEG 表
############################################################

DEG <- read.csv(
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/edgeR_TCX2-3-KO_vs_WT_sig.csv",
  header = TRUE
)


############################################################
### 5. 合并 DEG + RT 信息（先合一次，后面再补充 shift）
############################################################

DEG_RT <- DEG %>%
  left_join(gene_RT_full, by = c("X" = "Geneid"))


############################################################
### 6. 定义基因表达方向（Up / Down）
############################################################

DEG_RT <- DEG_RT %>%
  mutate(direction = case_when(
    logFC > 0 ~ "Up",
    logFC < 0 ~ "Down",
    TRUE      ~ "NS"
  ))


############################################################
### 7. 先把 RT 归一化：ESMS→MS，MSLS/ESLS/ESMSLS→LS，保留 Non-rep
############################################################

rt_collapse <- function(x){
  case_when(
    x %in% c("ES") ~ "ES",
    x %in% c("ESMS", "MS") ~ "MS",
    x %in% c("MSLS","LS","ESLS","ESMSLS") ~ "LS",
    x %in% c("Non-replication") ~ "Non-replication",
    TRUE ~ NA_character_
  )
}

gene_RT_full <- gene_RT_full %>%
  mutate(
    RT_WT_simple   = rt_collapse(RT_WT),
    RT_TCX2_simple = rt_collapse(RT_TCX2)
  )


############################################################
### 8. 在 ES/MS/LS 内部定义 Earlier / Later / No_change
###    对 Non-rep → S、S → Non-rep 单独分类
############################################################

# 只在 S 期内定义线性顺序
s_levels <- c("ES", "MS", "LS")

gene_RT_full <- gene_RT_full %>%
  mutate(
    WT_num   = ifelse(RT_WT_simple   %in% s_levels,
                      match(RT_WT_simple,   s_levels), NA_integer_),
    TCX2_num = ifelse(RT_TCX2_simple %in% s_levels,
                      match(RT_TCX2_simple, s_levels), NA_integer_)
  )

# 只在 S→S 的情况下计算数值 shift（用于散点图）
gene_RT_full <- gene_RT_full %>%
  mutate(
    RT_shift_num = ifelse(
      !is.na(WT_num) & !is.na(TCX2_num),
      TCX2_num - WT_num,
      NA_integer_
    )
  )

# 先定义 S 期内部的 Earlier / Later / No_change_S
gene_RT_full <- gene_RT_full %>%
  mutate(
    RT_shift_core = case_when(
      !is.na(WT_num) & !is.na(TCX2_num) & (TCX2_num < WT_num) ~ "Earlier",
      !is.na(WT_num) & !is.na(TCX2_num) & (TCX2_num > WT_num) ~ "Later",
      !is.na(WT_num) & !is.na(TCX2_num) & (TCX2_num == WT_num) ~ "No_change_S",
      TRUE ~ NA_character_
    )
  )

# 再考虑和 Non-replication 有关的情况
gene_RT_full <- gene_RT_full %>%
  mutate(
    RT_shift_type = case_when(
      # Non-rep → 进入 S 期（不称 earlier，单独记）
      RT_WT_simple == "Non-replication" & RT_TCX2_simple %in% s_levels ~ "NR_to_S",
      # S 期 → 变成 Non-replication
      RT_WT_simple %in% s_levels & RT_TCX2_simple == "Non-replication" ~ "S_to_NR",
      # 始终 Non-replication
      RT_WT_simple == "Non-replication" & RT_TCX2_simple == "Non-replication" ~ "NR_stable",
      # 其他 S→S 情况用核心分类
      !is.na(RT_shift_core) ~ RT_shift_core,
      TRUE ~ NA_character_
    )
  )

# 可选：设定 factor 顺序
gene_RT_full$RT_shift_type <- factor(
  gene_RT_full$RT_shift_type,
  levels = c("Earlier","Later","No_change_S","NR_to_S","S_to_NR","NR_stable")
)


############################################################
### 9. 把 RT shift 信息挂回 DEG_RT
############################################################

DEG_RT <- DEG_RT %>%
  left_join(
    gene_RT_full %>%
      dplyr::select(Geneid,
                    RT_WT_simple, RT_TCX2_simple,
                    RT_shift_num, RT_shift_type),
    by = c("X" = "Geneid")
  )


############################################################
### 10. 分类：RT shift × Expression（可选，用于着色）
############################################################

DEG_RT <- DEG_RT %>%
  mutate(category = case_when(
    RT_shift_type == "Earlier"    & direction == "Up"   ~ "Earlier & Up",
    RT_shift_type == "Earlier"    & direction == "Down" ~ "Earlier & Down",
    RT_shift_type == "Later"      & direction == "Up"   ~ "Later & Up",
    RT_shift_type == "Later"      & direction == "Down" ~ "Later & Down",
    RT_shift_type == "No_change_S" & direction == "Up"  ~ "No_change_S & Up",
    RT_shift_type == "No_change_S" & direction == "Down" ~ "No_change_S & Down",
    RT_shift_type %in% c("NR_to_S","S_to_NR","NR_stable") ~ as.character(RT_shift_type),
    TRUE ~ "Other"
  ))


############################################################
### 11. 可视化 1：RT shift_num × logFC（只看 S→S 的基因）
############################################################

ggplot(DEG_RT %>% filter(!is.na(RT_shift_num)),
       aes(x = RT_shift_num, y = logFC, color = category)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_classic() +
  xlab("RT shift within S-phase (TCX2_num - WT_num)") +
  ylab("logFC (TCX2 vs WT)") +
  ggtitle("RT shift (ES/MS/LS only) vs expression change")


############################################################
### 12. 可视化 2：按 RT_shift_type 看 logFC 分布（含 Non-rep 转变）
############################################################

ggplot(DEG_RT %>% filter(!is.na(RT_shift_type)),
       aes(x = RT_shift_type, y = logFC, fill = RT_shift_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~ direction, nrow = 2, scales = "free_y") +
  theme_classic() +
  ggtitle("Expression distribution across RT shift groups (Up vs Down)") +
  ylab("logFC") +
  xlab("RT shift type")


