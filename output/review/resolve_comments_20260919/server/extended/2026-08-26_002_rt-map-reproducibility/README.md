# 002 — Repli-ATAC RT 图谱重复性

## 目的

量化 `NIP.1`、`ZH11.1`、`ZH11.2` 三套用于最终 RT 分类的数据在 OCR 层面的重复性，并为最终类别提供区域级置信度。

## 输入

- `/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/replication_classification_results.csv`
- 当前目录只读软链接：`data/repli_atac_classification.csv`

## 方法

- 三套 phase call 两两计算 exact agreement、Cohen's κ 和阶段集合 Jaccard。
- 计算三套数据对最终分类的逐 OCR 支持度、是否全体一致、是否存在至少 2/3 完全一致。
- 本任务不重新分类，只评估现有分类的可重复性。

## 运行

```bash
python code/01_rt_map_reproducibility.py \
  --input data/repli_atac_classification.csv \
  --output-dir output --plot-dir plots
```

环境：Python 3.11.4；pandas 1.5.3；SciPy 1.10.1；scikit-learn 1.3.0；matplotlib 3.7.1。

## 主要结果

- 共评估 62,564 个 OCR；该数字与 PDF 正文的 62,468 相差 96，投稿前需核对版本。
- 三套数据完全一致的 OCR 占 51.26%；全部 OCR 都有至少 2/3 完全一致的 call。
- 两两 exact agreement 为 66.13%–68.68%，Cohen's κ 为 0.517–0.557。
- 最终分类与三套原始 call 的平均阶段集合 Jaccard 为 0.862。

结论：多数投票后的图谱具有可用的重复性，但约一半 OCR 并非三套数据完全一致。补充材料中应报告 region-level confidence，避免把所有 OCR 当成同等确定。

## 输出

- `output/rt_map_reproducibility_summary_002.tsv`
- `output/rt_map_pairwise_agreement_002.tsv`
- `output/rt_map_region_confidence_002.tsv.gz`
- `plots/rt_map_reproducibility_002.pdf`
- 日志：`output/logs/01_rt_map_reproducibility.log`

