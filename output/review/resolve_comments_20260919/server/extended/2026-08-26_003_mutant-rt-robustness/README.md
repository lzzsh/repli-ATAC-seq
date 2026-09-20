# 003 — CPP 突变体 RT 稳健性与等位基因一致性

## 目的

用未阈值化的 Repli-seq 原始 1-kb bin counts 检验 CPP8/CPP11 突变体的基线依赖 RT 重分布、独立等位基因一致性和 CUT&Tag-bound bins 的效应。

## 输入

- 原始 Repli-seq counts：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt`
- CPP8 bound bins：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/TCX2_overlap_segmentation.bed`
- CPP11 bound bins：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/SOL1_overlap_segmentation.bed`

## 方法

1. 每个 Repli-seq library 按 bin 长度和库深转为 TPM。
2. 同条件重复取均值；ES/MS/LS 分别除以对应 G1。
3. `WRT = (0.5 × MS + LS) / (ES + MS + LS)`；`ΔWRT = mutant − WT`。
4. `ΔWRT ~ WT WRT` 的正斜率定义为极化；1-Mb block bootstrap 和 cluster-robust SE。
5. 比较独立等位基因 ΔWRT 的 block-bootstrap Pearson correlation。

注意：早期测试曾错误使用为分类而截断的 `repli_phase_classified.csv` 信号，产生反向假象；该结果已被原始 counts 重跑结果完整覆盖。

## 运行

```bash
python code/01_mutant_rt_robustness.py \
  --input data/repli_raw_counts.txt \
  --cpp8-bins data/cpp8_bound_bins.bed \
  --cpp11-bins data/cpp11_bound_bins.bed \
  --output-dir output --plot-dir plots --bootstrap 1000
```

## 主要结果

- CPP8 两个等位基因均稳定极化：
  - `oscpp8-1` slope = 0.114，95% block-bootstrap CI 0.085–0.139。
  - `oscpp8-3` slope = 0.172，95% CI 0.141–0.198。
- CPP8 两等位基因 ΔWRT 高度一致：Pearson r = 0.858，95% CI 0.851–0.864；同号率 87.1%。
- CPP11 存在等位基因异质性：
  - `oscpp11-5` 为弱正极化，slope = 0.0426，95% CI 0.0220–0.0620。
  - `oscpp11-8` slope = −0.0515，95% CI −0.0723 至 −0.0311，不支持与 `oscpp11-5` 相同的极化方向。
  - 两等位基因的局部 ΔWRT 仍相关（r = 0.768），说明“总体图形相似”不等同于“基线依赖方向相同”。
- CUT&Tag-bound bins 的绝对 RT shift 略大于未结合背景，但没有显示更强的极化 score。

## 对稿件的直接影响

PDF 当前把 `oscpp8-3` 和 `oscpp11-8` 都描述为极化。服务器原始 counts、本任务结果和旧 `WT_vs_sol1_8_RT_quintile_boxplot_all_bins_stats.tsv` 均提示 `oscpp11-8` 不支持同等强度和方向的结论。正式投稿前必须核对样本身份、Figure 3J 数据源与作图脚本；在未解决前，机制主线应以 CPP8 为主，CPP11 写成等位基因依赖或较弱效应。

## 输出

- `output/mutant_polarization_effects_003.tsv`
- `output/independent_allele_concordance_003.tsv`
- `output/cuttag_bound_site_adjusted_effects_003.tsv`
- `plots/mutant_rt_robustness_003.pdf`

