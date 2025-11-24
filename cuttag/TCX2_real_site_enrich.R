setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/cuttag_tcx/macs2/macs2_p1e-5/tcx_cut_out/")
tcx2_site <- read.table("./TCX2_overlap_peak.bed")

peaks_reads <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_org.gff3", sep = "\t")
colnames(peaks_reads) <- c("chr","start","end","RT")

# Genomic RT length distribution
RT_freq_genome <- peaks_reads %>% 
  group_by(RT) %>% 
  summarise(Sum = sum(end) - sum(start)) %>%
  filter(!(RT %in% c("ESLS","ESMSLS"))) %>%
  mutate(
    percent = Sum / sum(Sum),
    RT = factor(RT, levels = c("ES","ESMS","MS","MSLS","LS"))
  )

# Observed TCX2 peaks per RT
tcx2_RT_table <- RT_freq_tcx2 %>% dplyr::select(RT, Count)

# Expected counts from genomic proportions
total_tcx2 <- sum(tcx2_RT_table$Count)
RT_expected <- RT_freq_genome %>%
  dplyr::select(RT, percent) %>%
  mutate(Expected = percent * total_tcx2)

# Combine observed vs expected
RT_compare <- RT_expected %>%
  left_join(tcx2_RT_table, by="RT") %>%
  mutate(Observed = Count) %>%
  dplyr::select(RT, Observed, Expected)

print(RT_compare)

# Chi-square test for global deviation
chisq_test <- chisq.test(RT_compare$Observed, p = RT_freq_genome$percent)
chisq_test

# Fisher test per RT class
library(purrr)

fisher_results <- map_df(RT_freq_genome$RT, function(rt_class) {
  
  observed <- RT_compare$Observed[RT_compare$RT == rt_class]
  expected <- RT_compare$Expected[RT_compare$RT == rt_class]
  
  observed_other <- sum(RT_compare$Observed) - observed
  expected_other <- sum(RT_compare$Expected) - expected
  
  m <- matrix(c(observed, observed_other,
                expected, expected_other), nrow=2, byrow=TRUE)
  
  data.frame(
    RT = rt_class,
    observed = observed,
    expected = expected,
    fisher_p = fisher.test(m)$p.value,
    enrichment = observed / expected
  )
})

print(fisher_results)