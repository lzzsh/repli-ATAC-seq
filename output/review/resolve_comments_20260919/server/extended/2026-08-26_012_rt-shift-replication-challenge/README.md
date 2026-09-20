# 012 — 两个CPP8等位基因共同RT异常与复制挑战特征

## 状态

**正式分析已完成（2026-08-26）。** 目录使用`replication-challenge`，因为没有位点分辨率DNA损伤或复制叉测量。

## 科学问题与结论

检验两个CPP8等位基因共同RT异常是否位于具有潜在复制挑战的基因组环境，并与更大的ATAC/RNA扰动共定位。

匹配分析显示，top-5%共识RT异常bins富集Class I TE，靠近centromere-repeat proxy，并位于WT表达、gene density和OCR density较低的环境；Class II TE和最近基因长度在primary匹配中不稳定。热点OCR的绝对ΔATAC在两alleles三阶段均增大，MS最强；邻近基因的绝对RNA logFC也小幅增加。

## 输入与注释

- RT bins：`../2026-08-26_009_cpp8-allele-path-replication/output/repli_rt_bin_metrics_009.tsv.gz`
- raw Repli-seq：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt`
- OCR metrics：`../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz`
- gene metrics：`../2026-08-26_009_cpp8-allele-path-replication/output/gene_level_path_metrics_009.tsv.gz`
- RNA joint model：`../2026-08-26_009_cpp8-allele-path-replication/output/rna_joint_edger_all_genes_009.tsv.gz`
- FASTA：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta`
- GFF3：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3`
- TE annotation：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/annotate.bed`；Class I 102,687条，Class II 291,003条。
- Histone BED：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/histone/`中的H3K27ac、H3K4me1、H3K4me3、H3K9me2和H3K27me3。

`centromere-repeat proxy`由每条染色体GFF3中含`centromere-specific`注释的基因/retrotransposon cluster midpoint中位数定义，不是官方centromere/pericentromere边界。

## 共识RT异常定义

```text
若两个alleles的ΔRT同号：
  consensus_shift = sign(mean ΔRT) × min(|ΔRT1|, |ΔRT3|)
否则：
  consensus_shift = 0
consensus_instability = |consensus_shift|
```

- 373,251个1-kb bins进入构建。
- top-5% threshold=`0.2632919444`；类别为stable 295,115、discordant 59,473、concordant later 17,966、concordant earlier 697。
- 阈值敏感性：top 2.5%、5%、10%。

## 统计设计

连续多变量模型调整WT WRT spline、WT G1 coverage和染色体，并同时评估8个primary features。热点匹配在以下strata内为每个hotspot抽取最多5个control：

- chromosome；
- WT WRT decile；
- GC quintile；
- G1 coverage quintile；
- gene/TE compartment。

匹配差异以1-Mb blocks进行1,000次bootstrap；seed=`20260826`。每个阈值内8个features进行BH FDR。

## 正式结果

### Top-5% matched hotspot analysis

| feature | standardized hotspot−control | 95% block-bootstrap CI | FDR |
|---|---:|---:|---:|
| Class I TE fraction | 0.1507 | 0.1252–0.1763 | 0.00266 |
| distance to centromere-repeat proxy | −0.4989 | −0.6187至−0.3838 | 0.00266 |
| WT expression within 50 kb | −0.5177 | −0.5862至−0.4494 | 0.00266 |
| gene density within 50 kb | −0.1277 | −0.1886至−0.0713 | 0.00266 |
| OCR density within 50 kb | −0.2549 | −0.3154至−0.2059 | 0.00266 |
| GC fraction | −0.0082 | −0.0126至−0.0039 | 0.00266 |
| Class II TE fraction | 0.0129 | −0.0135–0.0403 | 0.344 |
| nearest-gene length | −0.0201 | −0.0433–0.0034 | 0.107 |

Class I TE、centromere-proxy距离、低WT表达、低gene/OCR density在2.5%、5%和10%阈值方向稳定。最近基因长度仅在10%阈值出现很弱负差异，不能支持“long-gene hotspot”结论。

### 连续模型

在344,742个完整bins中，Class I TE β=0.0812，centromere-proxy距离β=−0.1909，WT表达β=−0.1539；方向与匹配分析一致。Class II TE在连续模型中为小正效应，但matched primary结果不显著，因此不作为headline。

### ATAC与RNA连接

| layer | allele/stage | standardized hotspot difference | 95% cluster CI |
|---|---|---:|---:|
| `|ΔATAC|` | `oscpp8-1` MS | 0.1323 | 0.0859–0.1787 |
| `|ΔATAC|` | `oscpp8-3` MS | 0.2083 | 0.1672–0.2494 |
| `|RNA logFC|` | `oscpp8-1` MS, ±1 kb genes | 0.0730 | 0.0085–0.1375 |
| `|RNA logFC|` | `oscpp8-3` MS, ±1 kb genes | 0.0912 | 0.0223–0.1601 |

两alleles的ES/MS/LS OCR `|ΔATAC|`差异均为正，MS最强。

## 运行顺序

```bash
python3 code/01_build_consensus_and_features.py \
  --rt-bins ../2026-08-26_009_cpp8-allele-path-replication/output/repli_rt_bin_metrics_009.tsv.gz \
  --raw-repli /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt \
  --fasta /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta \
  --gff /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3 \
  --annotation-bed /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/annotate.bed \
  --histone-dir /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/histone \
  --ocr-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz \
  --rna ../2026-08-26_009_cpp8-allele-path-replication/output/rna_joint_edger_all_genes_009.tsv.gz \
  --output-dir output
python3 code/03_fit_feature_models.py \
  --features output/replication_challenge_features_012.tsv.gz --output-dir output
python3 code/04_hotspot_matched_and_links.py \
  --features output/replication_challenge_features_012.tsv.gz \
  --ocr-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz \
  --gene-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/gene_level_path_metrics_009.tsv.gz \
  --output-dir output --bootstrap 1000
python3 code/05_plot_replication_challenge.py \
  --models output/feature_multivariable_models_012.tsv \
  --matched output/matched_hotspot_enrichment_012.tsv \
  --signed output/signed_shift_class_models_012.tsv --plot-dir plots
```

## 软件环境

Python 3.11.4；NumPy 1.26.4；pandas 1.5.3；SciPy 1.10.1；statsmodels 0.14.0；matplotlib 3.7.1；pysam 0.22.1。

## 稿件位置与解释边界

- 若Figure 8已有γH2AX/DNA-repair主线，可加入一个matched-feature panel作为“global damage response”与“RT-sensitive genomic context”的桥接。
- 完整阈值、histone、signed-class和多变量模型进入SI。
- 可写`replication-challenge-associated genomic context`；不可写`DNA-damage hotspots`、`common fragile sites`、`fork stalling sites`或“CPP8直接保护这些位点”。

## 输出

- 表：`cpp8_consensus_rt_instability_012.tsv.gz`、`replication_challenge_features_012.tsv.gz`、`feature_multivariable_models_012.tsv`、`matched_hotspot_enrichment_012.tsv`、`signed_shift_class_models_012.tsv`、`hotspot_atac_rna_links_012.tsv`。
- 图：`rt_instability_feature_forest_012.pdf`、`hotspot_matched_enrichment_012.pdf`、`signed_shift_feature_profiles_012.pdf`。
