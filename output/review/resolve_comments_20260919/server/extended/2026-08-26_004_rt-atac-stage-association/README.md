# 004 — OCR 层面 ΔRT–分阶段 ΔATAC 关联

## 目的

检验遗传扰动引起的 RT shift 是否与 ES、MS、LS 阶段的局部可及性变化相关，并控制空间相关和基线状态。

## 输入

- ATAC OCR counts：`/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/signal/ATAC_peaks_quan_peak.txt`
- Repli-seq raw counts：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt`

ATAC 列顺序由 `ATAC-seq-CR-2/signal/signal_compare.sh` 确认。Repli-seq 的样本顺序由 `segmentation/01_quan.sh` 确认。

## 方法

- ATAC counts 使用 positive-count median-ratio size factors；ES 重复取归一化均值，MS/LS 为单样本。
- Repli-seq 使用 TPM、重复均值、S/G1 和 WRT 公式。
- OCR 按真实 overlap bp 加权映射到 1-kb RT bins。
- 每个突变体、阶段分别拟合：`ΔATAC ~ ΔRT + WT WRT + WT ATAC + OCR length + chromosome`。
- 推断采用 1-Mb block cluster-robust SE 和 1,000 次 block bootstrap。

## 主要结果

- 89,361 个 OCR 中，约 88.8k 个获得可用的连续 RT 映射。
- CPP8：ΔRT–ΔATAC 在三个阶段均为正，MS 最强。
  - ES β = 0.0484，95% bootstrap CI 0.0358–0.0603。
  - MS β = 0.1428，95% CI 0.1308–0.1536。
  - LS β = 0.1099，95% CI 0.0955–0.1218。
- CPP8 并非 ES→MS→LS 单调增强，而是在 MS 达峰；线性阶段趋势不显著（P = 0.285）。
- CPP11：ES 不显著；MS β = 0.0241、LS β = 0.0297，效应较小但通过空间稳健检验。

结论：CPP8 RT redistribution 与同一 OCR 的阶段化可及性改变具有稳定但非单调的关联，最明显发生在 MS；这比“全局可及性重塑”更符合局部、阶段依赖的模型。

## 输出

- `output/ocr_rt_atac_metrics_004.tsv.gz`
- `output/rt_atac_stage_associations_004.tsv`
- `output/rt_atac_es_to_ls_trend_tests_004.tsv`
- `plots/rt_atac_stage_association_004.pdf`

