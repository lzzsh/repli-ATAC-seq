# 009 — CPP8两个等位基因的完整RT–ATAC–RNA路径复现

## 状态

**正式分析已完成（2026-08-26）。** 本目录仅包含`oscpp8-1`与`oscpp8-3`的双等位基因复现，不混入CPP11、阶段匹配、空间尺度或复制挑战分析。

## 科学问题与结论

问题是：CPP8缺失相关的`ΔRT → promoter ΔATAC → RNA logFC`统计路径能否在两个独立等位基因中重复。

结论是：两个等位基因的RT重分布高度一致，OCR层面的ΔRT–ΔATAC关联在ES、MS和LS均为正；MS启动子路径在±1 kb和±3 kb均重复。路径方向一致，但`oscpp8-1`的间接效应显著小于`oscpp8-3`，因此应写成“genetic replication with allele-dependent effect size”，而不是完全同质的机制复现。

## 输入

- Repli-seq raw counts：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt`
- ATAC sample sheet：`/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/sampleSheet1.csv`
- 固定OCR集合：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/peaks_quan_tpm_filtered_pos_org.txt`
- RNA raw counts：`/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/gene_counts.txt`
- FASTA：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta`
- GFF3：`/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3`

## 数据重建

1. 用`bedtools multicov`在89,361个固定OCR上重新计数WT、TCX2_1和TCX2_3共12个ATAC BAM；12个library联合计算positive-count median-ratio size factors。
2. RNA从55,986个raw-count genes重新建模。分别对WT–`oscpp8-1`和WT–`oscpp8-3`执行`filterByExpr`，取交集后在六样本联合edgeR QL模型中估计两个contrast，得到24,212个共同可检测基因。
3. Repli-seq由raw counts计算TPM、重复均值、S/G1与WRT；`ΔRT=WRTmutant−WRTWT`。
4. GFF3的gene条目和链方向用于重建TSS；±1 kb映射22,043个路径基因，±3 kb映射23,460个路径基因。

旧`edgeR_TCX2-*-KO_vs_WT.csv`实际仅保留FDR约0.05的显著基因，不能作为无结局筛选的全基因路径集合。旧分析的`a×b≈0.033`因此不再作为正式主结果；本任务基于共同可检测基因的重算值是稿件应采用的版本。

## 模型与推断

- OCR模型：`z(ΔATAC) ~ z(ΔRT) + WT WRT + WT stage ATAC + log(OCR length+1) + chromosome`。
- promoter a路径：`z(ΔATAC) ~ z(ΔRT) + pretreatment covariates + chromosome`。
- promoter b/c′路径：`z(RNA logFC) ~ z(ΔATAC) + z(ΔRT) + pretreatment covariates + chromosome`。
- pretreatment covariates：WT WRT、WT ATAC、WT RNA abundance、启动子OCR数量。
- 推断单位：1-Mb genomic blocks；379个有效blocks。
- OCR bootstrap：1,000次；路径同步bootstrap：2,000次；空间null：500次。
- promoter窗口：TSS ±1 kb和±3 kb；随机种子：`20260826`。
- 12个间接效应使用BH FDR；allele-difference contrasts单独校正。

## 正式结果

### 等位基因一致性

- 360,293个共同1-kb bins的ΔRT：Pearson `r=0.857929`，95% block-bootstrap CI 0.851228–0.864614；Spearman `ρ=0.874680`。
- RNA logFC跨allele相关：`r=0.628356`，95% CI 0.609877–0.646061。
- promoter ΔATAC跨allele相关为0.549–0.663，说明下游方向可重复但效应并非完全一致。

### OCR层面ΔRT–ΔATAC

| allele | ES β (95% CI) | MS β (95% CI) | LS β (95% CI) |
|---|---:|---:|---:|
| `oscpp8-1` | 0.0538 (0.0417–0.0655) | 0.1262 (0.1150–0.1367) | 0.1216 (0.1079–0.1339) |
| `oscpp8-3` | 0.0481 (0.0349–0.0608) | 0.1432 (0.1310–0.1546) | 0.1094 (0.0957–0.1229) |

六项区块bootstrap CI均排除0；两个等位基因都支持局部RT–accessibility coupling。

### MS启动子路径

| allele | window | a | b | a×b | 95% block-bootstrap CI | FDR |
|---|---:|---:|---:|---:|---:|---:|
| `oscpp8-1` | ±1 kb | 0.1053 | 0.1551 | 0.01634 | 0.01249–0.02019 | 0.00120 |
| `oscpp8-3` | ±1 kb | 0.1191 | 0.2231 | 0.02657 | 0.02229–0.03085 | 0.00120 |
| `oscpp8-1` | ±3 kb | 0.1172 | 0.1569 | 0.01839 | 0.01454–0.02233 | 0.00120 |
| `oscpp8-3` | ±3 kb | 0.1257 | 0.2204 | 0.02771 | 0.02301–0.03223 | 0.00120 |

LS也在两alleles和两个窗口中保持正向且CI排除0；ES中`oscpp8-1`较弱，±3 kb CI跨0。

### 效应异质性

- MS ±1 kb：`oscpp8-1 − oscpp8-3 = −0.01023`，95% CI −0.01474至−0.00571。
- MS ±3 kb：差值`−0.00933`，95% CI −0.01406至−0.00451。

因此“方向和路径重复”成立，但效应大小具有allele dependence。

### LOCO与空间null

- MS ±1 kb，加入ΔATAC相对baseline的pooled LOCO ΔR²：`oscpp8-1=0.02165`，`oscpp8-3=0.03653`。
- 加入ΔRT的对应ΔR²：`oscpp8-1=0.00056`，`oscpp8-3=−0.00007`。
- MS两个alleles在500次空间null中的经验`P=0.001996`；±1 kb与±3 kb结论一致。

## 运行顺序

在本任务目录执行：

```bash
bash code/01_build_joint_atac_matrix.sh
Rscript code/02_build_joint_rna_contrasts.R
python3 code/03_cpp8_allele_path_replication.py \
  --atac-counts output/atac_all_cpp8_alleles_counts_009.tsv.gz \
  --repli /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt \
  --gff /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3 \
  --rna output/rna_joint_edger_all_genes_009.tsv.gz \
  --output-dir output --plot-dir plots \
  --bootstrap-ocr 1000 --bootstrap-path 2000 --windows 1000,3000
python3 code/04_cpp8_allele_validation.py \
  --gene-metrics output/gene_level_path_metrics_009.tsv.gz \
  --output-dir output --plot-dir plots --permutations 500
```

## 软件环境

- Python 3.11.4；NumPy 1.26.4；pandas 1.5.3；SciPy 1.10.1；statsmodels 0.14.0；scikit-learn 1.3.0；matplotlib 3.7.1。
- R 4.4.3；edgeR 4.2.1。
- bedtools 2.30.0。

## 解释边界

- 三种组学来自独立bulk样本集合，不能称为严格causal mediation。
- 1-Mb block inference控制空间伪重复，但不能替代生物学重复。
- RT对RNA的独立全局预测贡献很小；证据支持局部多组学耦联，不支持“RT在全基因组直接驱动转录”。
- 不报告mediated proportion，因为总效应接近0且存在competitive/suppression-like structure。

## 关键输出

- 数据层：`atac_all_cpp8_alleles_counts_009.tsv.gz`、`rna_joint_edger_all_genes_009.tsv.gz`、`repli_rt_bin_metrics_009.tsv.gz`、`ocr_rt_atac_metrics_009.tsv.gz`、`gene_level_path_metrics_009.tsv.gz`。
- 统计表：`cpp8_allele_ocr_associations_009.tsv`、`cpp8_allele_path_effects_009.tsv`、`cpp8_allele_path_heterogeneity_009.tsv`、`cpp8_allele_concordance_009.tsv`、`cpp8_allele_loco_009.tsv`、`cpp8_allele_spatial_null_009.tsv`。
- 图：`plots/cpp8_allele_path_replication_009.pdf`、`plots/cpp8_allele_validation_009.pdf`。
- `*_smoke.log`仅为低迭代调试记录，正式数值以无`_smoke`日志和当前TSV为准。
