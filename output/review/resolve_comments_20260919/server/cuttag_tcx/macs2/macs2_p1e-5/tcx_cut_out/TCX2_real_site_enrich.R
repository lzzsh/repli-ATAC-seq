# =========================
# TCX2 peaks ~ RT enrichment
# Fisher background: number of RT segments (peaks_reads$RT counts), NOT genome bp
# Background from: table(peaks_reads$RT)
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# -------- paths --------
setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/cuttag_tcx/macs2/macs2_p1e-5/tcx_cut_out/")

tcx2_bed <- "./TCX2_overlap_peak.bed"
rt_gff3  <- "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_org.gff3"

# -------- read --------
tcx2_site <- read.table(tcx2_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

peaks_reads <- read.table(rt_gff3, sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(peaks_reads)[1:4] <- c("chr","start","end","RT")
peaks_reads <- peaks_reads[, c("chr","start","end","RT")]

# -------- settings --------
RT_keep <- c("ES","ESMS","MS","MSLS","LS")  # 你要分析的5类

# -------------------------
# 1) Background: RT segment counts (table(peaks_reads$RT))
# -------------------------
RT_freq_bg <- as.data.frame(table(peaks_reads$RT), stringsAsFactors = FALSE) %>%
  rename(RT = Var1, BgCount = Freq) %>%
  mutate(RT = as.character(RT)) %>%
  filter(RT %in% RT_keep) %>%
  mutate(RT = factor(RT, levels = RT_keep)) %>%
  arrange(RT) %>%
  mutate(
    percent_bg = BgCount / sum(BgCount)   # 背景概率：各RT段数占比
  )

print(RT_freq_bg)

# -------------------------
# 2) Observed: TCX2 overlap counts per RT (table(tcx2_site$V4))
# -------------------------
RT_freq_tcx2 <- as.data.frame(table(tcx2_site$V4), stringsAsFactors = FALSE) %>%
  rename(RT = Var1, Count = Freq) %>%
  mutate(RT = as.character(RT)) %>%
  filter(RT %in% RT_keep) %>%
  mutate(RT = factor(RT, levels = RT_keep)) %>%
  arrange(RT)

print(RT_freq_tcx2)

# -------------------------
# 3) Expected counts from background segment proportions
# -------------------------
total_tcx2 <- sum(RT_freq_tcx2$Count)

RT_expected <- RT_freq_bg %>%
  select(RT, BgCount, percent_bg) %>%
  mutate(Expected = percent_bg * total_tcx2)

# -------------------------
# 4) Combine observed vs expected
# -------------------------
RT_compare <- RT_expected %>%
  left_join(RT_freq_tcx2, by = "RT") %>%
  mutate(Observed = Count) %>%
  select(RT, Observed, Expected, BgCount, percent_bg) %>%
  arrange(RT)

print(RT_compare)

# -------------------------
# 5) Chi-square test (global deviation) using background segment proportions
# -------------------------
chisq_test <- chisq.test(RT_compare$Observed, p = RT_compare$percent_bg)
print(chisq_test)

# -------------------------
# 6) Fisher test per RT class using background segment counts
#    2×2 table for each RT:
#      [ TCX2_in_RT   TCX2_other ]
#      [ BG_in_RT     BG_other   ]
# -------------------------
fisher_results <- map_df(levels(RT_compare$RT), function(rt_class) {
  
  a <- RT_compare$Observed[RT_compare$RT == rt_class]
  b <- sum(RT_compare$Observed) - a
  
  c <- RT_compare$BgCount[RT_compare$RT == rt_class]
  d <- sum(RT_compare$BgCount) - c
  
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  colnames(m) <- c("This_RT", "Other_RT")
  rownames(m) <- c("TCX2_peaks", "RT_segments_bg")
  
  ft <- fisher.test(m)
  
  # enrichment: (TCX2比例)/(背景比例)
  enrich <- (a / sum(RT_compare$Observed)) / (c / sum(RT_compare$BgCount))
  
  data.frame(
    RT = as.character(rt_class),
    observed = as.numeric(a),
    expected = RT_compare$Expected[RT_compare$RT == rt_class],
    fisher_p = ft$p.value,
    enrichment = as.numeric(enrich),
    odds_ratio = unname(ft$estimate),
    stringsAsFactors = FALSE
  )
})

print(fisher_results)
print(fisher_results %>% arrange(fisher_p))

# -------------------------
# 7) Contingency tables for Supplementary (one row per RT)
# -------------------------
contingency_table <- map_df(levels(RT_compare$RT), function(rt_class) {
  
  a <- RT_compare$Observed[RT_compare$RT == rt_class]
  b <- sum(RT_compare$Observed) - a
  
  c <- RT_compare$BgCount[RT_compare$RT == rt_class]
  d <- sum(RT_compare$BgCount) - c
  
  data.frame(
    RT = as.character(rt_class),
    TCX2_in_RT = as.numeric(a),
    TCX2_other = as.numeric(b),
    BG_in_RT = as.numeric(c),
    BG_other = as.numeric(d),
    stringsAsFactors = FALSE
  )
})

print(contingency_table)
