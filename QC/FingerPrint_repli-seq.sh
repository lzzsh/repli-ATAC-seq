# FingerPrint
plotFingerprint \
  -b bwa_all_rawdata/Rep-T3-1_S99_1_sort_rmdup_rmor_q30.bam bwa_all_rawdata/Rep-T3-2_S99_1_sort_rmdup_rmor_q30.bam bwa_all_rawdata/Rep-T3-3_S99_1_sort_rmdup_rmor_q30.bam \
  bwa_all_rawdata/Rep-T3-4_S99_1_sort_rmdup_rmor_q30.bam bwa_all_rawdata/Rep-T3-5_S99_1_sort_rmdup_rmor_q30.bam \
  /storage/liuxiaodongLab/liaozizhuo/Projects/repli-seq-pre/bwa_all_rawdata/Rep-T2-2_S99_1_sort_rmdup_rmor_q30.bam \
  /storage/liuxiaodongLab/liaozizhuo/Projects/repli-seq-pre/bwa_all_rawdata/Rep-T2-2_S99_1_sort_rmdup_rmor_q30.bam \
  /storage/liuxiaodongLab/liaozizhuo/Projects/repli-seq-pre/bwa_all_rawdata/Rep-T2-2_S99_1_sort_rmdup_rmor_q30.bam \
  -o qc_plots/fingerprint_signal_distribution_all.pdf \
  --plotTitle "Signal Distribution Fingerprint"