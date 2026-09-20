# 006 — 因果路径的空间验证、置换与跨染色体预测

## 目的

独立验证任务 005 的路径是否具有跨染色体泛化能力，是否对启动子窗口稳健，以及是否超出保留染色体/局部空间结构的负对照。

## 输入

- 任务 005：`output/gene_level_path_metrics_005.tsv.gz`
- 任务 005：`output/rt_atac_rna_path_effects_005.tsv`

## 方法

### Chromosome leave-one-out prediction

每次留出一条染色体，在其余 11 条训练普通线性模型，并在未见过的染色体预测 RNA logFC：

- baseline：WT WRT、WT ATAC、WT RNA、OCR 数。
- baseline + ΔRT。
- baseline + ΔATAC。
- baseline + ΔRT + ΔATAC。

以 pooled out-of-chromosome R² 和 12 条染色体的 paired Wilcoxon test 比较模型。

### 负对照

每个突变体 × 阶段 × promoter window 做 500 次：

- 染色体内、WT WRT/WT ATAC 五分位分层置换。
- 每条染色体内按 TSS 排序后的 circular shift，保留局部自相关结构但打断 locus 对齐。

## 主要结果

- ΔATAC 在所有 12 个组合中均显著提升跨染色体 RNA 预测：ΔR² = 0.0367–0.0748；染色体配对检验均支持改善。
- ΔRT 单独的 ΔR² 基本为 0，范围约 −0.00029 至 +0.00002。
- 已知 ΔATAC 后再加入 ΔRT，增益仍接近 0；CPP8 MS 的数值增益最大但仅约 0.0019–0.0020，染色体层面不稳定。
- CPP8 MS 和 LS 的间接效应在 ±1/±3 kb 均保持同号、bootstrap CI 排除 0，并在四种空间负对照下 empirical P = 0.001996，scheme-specific FDR < 0.005。
- CPP8 ES 虽两个窗口的原始 CI 均排除 0，但多重校正与窗口一致性较边缘，不建议作为主结论。

## 解释

跨染色体预测明确把 ΔATAC 定位为与 RNA 改变最直接、最可泛化的组学层。ΔRT 对 RNA 的全局独立预测贡献不足，但其与 CPP8 MS/LS promoter ΔATAC 的空间共定位通过了严格负对照。最稳妥的结论是局部路径证据，而非全局 RT→RNA 因果效应。

## 输出

- `output/chromosome_loco_cv_summary_006.tsv`
- `output/chromosome_loco_cv_model_comparisons_006.tsv`
- `output/causal_path_permutation_tests_006.tsv`
- `output/promoter_window_sensitivity_006.tsv`
- `plots/causal_path_validation_006.pdf`
- 正式运行日志：`output/logs/01_causal_path_validation.log`（500 permutations；4:30.66；峰值内存约 450 MB）

