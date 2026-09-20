# 015 — Manuscript figure and analysis audit

## 目的

逐图检查 `Repli-ATAC-seq-20260705.pdf` 中 Figure 1–8 和 Figure S1–S5 的证据链、统计设计、重复结构、图注、Methods、源数据与可复现性，并对高风险数值进行独立复算。

本任务采用 `nature-writing` 的 Results/Discussion 审计原则：区分 observation、association、mechanism 与 causality；每项主张均检查 claim、figure、caption、source、script、statistics、replicates、threshold 和 genome version。

## 输入

- 稿件：`/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/Repli-ATAC-seq-20260705.pdf`
- FASTA：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta`
- GFF3：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3`
- 其余输入均来自项目根目录下既有只读分析目录，完整路径记录在脚本、逐图审计表和综合报告中。

输入文件 SHA-256：

- manuscript PDF：`2e704cd7c68e0d19000c497c9b826d72a23e5c54ab060d36b1356f73c17d0750`
- FASTA：`751b4091ad78d4ac0fd64432573bd91af0d286277f4e52b6ba73fd2acd12aecd`
- GFF3：`6dcc949fe878f8218ff2d8afd4bc0564c8560acbffa36e62b90b1e323489491a`

## 稿件页码

- Main figures：Figure 1 page 4；2 page 7；3 page 11；4 page 15；5 page 18；6 page 20；7 page 22；8 page 25。
- Supplementary figures：S1 page 36；S2 page 38；S3 page 40；S4 page 41；S5 page 43。

## 执行

PDF 页面和文本已输出至 `output/pdf_pages/` 与 `output/page_text/`。高风险面板复算：

```bash
python3 code/02_targeted_figure_reproduction.py \
  > output/logs/02_targeted_figure_reproduction.log 2>&1
```

环境：Python 3.11.4；pandas 1.5.3；SciPy 1.10.1。随机重采样类独立分析沿用已有任务的 seed `20260826` 和 1-Mb genomic blocks。

## 主要输出

- `output/MANUSCRIPT_FIGURE_AUDIT_015.md`：完整逐图结论和整改优先级。
- `output/figure_by_figure_audit_015.tsv`：Figure 1–8、S1–S5 的审计矩阵。
- `output/critical_submission_blockers_015.md`：投稿前阻断项。
- `output/figure1_ocr_count_check_015.tsv`：Figure 1 OCR 数量和比例。
- `output/figure2_variant_density_check_015.tsv`、`figureS2_maf_spectrum_check_015.tsv`：Figure 2/S2 GWAS复核。
- `output/figure3_polarization_check_015.tsv`：Figure 3 raw-count WRT方向复核。
- `output/figure7_dar_counts_015.tsv`、`figure7_peak_overlap_check_015.tsv`：Figure 7 DAR和Venn复核。
- `output/figure8_atac_rna_regression_check_015.tsv`、`figure8_overlap_universe_sensitivity_015.tsv`：Figure 8 R²与overlap复核。

## 可复现性分级

- 数值可复现：Figure 1 class composition、Figure 2 variant density/MAF、Figure 7 DAR counts/Venn、Figure 8 R²/Venn。
- 可从 raw counts 独立重建且发现矛盾：Figure 3 WRT polarization。
- 部分可复现：Figure 5 Wox11 RT、Figure 6 selected RT tracks、Figure 7/8多组学部分。
- 当前不可复现：Figure 4完整模型；Figure 5 WOX11 RNA/GO；Figure 6 phenotype；Figure 8 blot。

## 解释边界

- bulk RT、ATAC 和 RNA 来自不同样本集合，不能写成严格 causal mediation。
- 1-Mb genomic-block inference控制空间伪重复，但不能替代缺失的 biological replicates。
- CUT&Tag来自 ZS97 leaf-sheath protoplasts，RT/ATAC/RNA主要来自 ZH11 root tips；结合位点应写成 cross-context association，不能自动等同 root-tip direct targets。
