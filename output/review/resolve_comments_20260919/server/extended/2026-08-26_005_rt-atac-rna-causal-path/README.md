# 005 — 启动子层面 RT→ATAC→RNA 路径分析

## 目的与表述边界

估计由 CPP8/CPP11 遗传扰动锚定的 locus-level `ΔRT → promoter ΔATAC → RNA logFC` 路径。由于三种组学来自不同样本集合，本任务提供的是 causal-path-compatible evidence，不是同一个体内的严格中介因果证明。

## 参考基因组与注释

- 已确认参考基因组：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta`
- 已确认注释：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3`
- TSS 从 GFF3 的 `gene` 条目、链方向和 `ID=LOC_Os...` 重新生成。未使用第 4 列缺少基因 ID 的旧 `TSS.bed`。

## 其他输入

- OCR ATAC metrics：任务 004 的 `output/ocr_rt_atac_metrics_004.tsv.gz`
- gene-level RT：`RNA-seq/star_all_rawdata/featureCounts_out/gene_TSS_mapped_to_repli_bin.unique_gene_level.txt`
- RNA edgeR：`edgeR_TCX2-3-KO_vs_WT.csv`、`edgeR_SOL1-8-KO_vs_WT.csv`

gene-level RT 文件由旧 `WRT_RNA-seq.R` 从未阈值化 Repli-seq counts 计算；RNA outcome 使用完整 edgeR logFC，而不是只保留显著 DEG。

## 方法

- 启动子敏感性窗口：TSS ±1 kb、±3 kb。
- 同一启动子内 OCR 的归一化 counts 求和后计算 ATAC log2 fold change。
- a 路径：`ΔATAC ~ ΔRT + pretreatment covariates + chromosome`。
- b/c′ 路径：`RNA logFC ~ ΔATAC + ΔRT + pretreatment covariates + chromosome`。
- pretreatment covariates：WT WRT、WT ATAC、WT RNA abundance、启动子 OCR 数量。
- 1-Mb block cluster-robust SE；2,000 次同步 block bootstrap；间接效应为标准化 `a × b`。
- 不报告在总效应接近 0 时不稳定的 mediated proportion。

## 主要结果

### CPP8

- MS 是最稳定路径：
  - ±1 kb：a×b = 0.0331，95% CI 0.0239–0.0425，FDR = 0.0060。
  - ±3 kb：a×b = 0.0330，95% CI 0.0237–0.0431，FDR = 0.0060。
- LS 路径较弱但跨窗口稳定：
  - ±1 kb：0.0112，95% CI 0.00358–0.0188，FDR = 0.0120。
  - ±3 kb：0.0130，95% CI 0.00393–0.0235，FDR = 0.0080。
- ES 为边缘效应；±1 kb FDR 约 0.050，±3 kb FDR 约 0.072。
- b 路径在所有阶段都强且为正：promoter ΔATAC 与 RNA logFC 的标准化系数约 0.22–0.30。

### CPP11

- 没有跨 ±1/±3 kb 稳定、经多重校正的正向间接路径。
- ±3 kb ES 出现弱负间接效应（FDR 约 0.050），与 `oscpp11-8` 缺乏稳定极化的 QC 结果一致，不应与 CPP8 合并解释。

### 总效应与竞争路径

- CPP8 的 RT→RNA 总效应大多接近 0且不显著。
- 在 MS 中，正向间接效应与负向 c′ 同时存在，属于 competitive/suppression-like structure。
- 因此不能写成“RT shift 在全基因组直接驱动 RNA 改变”；更准确的是“CPP8 perturbation 下，RT shift 与 promoter accessibility change 共定位，而 accessibility 与 transcriptional output 紧密相连”。

## 输出

- `output/rt_atac_rna_path_effects_005.tsv`
- `output/gene_level_path_metrics_005.tsv.gz`
- `output/gene_tss_from_all_diy_gff3_005.tsv`
- `plots/rt_atac_rna_causal_paths_005.pdf`

