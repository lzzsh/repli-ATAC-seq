# 014 — DNA-damage/cell-cycle transcriptional response locality

## 目的

检验论文 Figure 8 中 DNA repair、SOG1 target 和 cell-cycle transcriptional response 是否优先定位于 CPP8 缺失后发生较强局部 RT/ATAC 扰动的基因，而不是仅表现为全局 stress-response programme。

本任务是基因集层面的局部性检验，不是严格的个体水平因果中介分析。

## 预定义基因集

- `manuscript_dna_repair`：论文及旧分析使用的 9 个 DNA-repair panel genes。该集合由 RNA 结果选择，因此 RNA logFC 仅作描述性正对照；正式检验仅评价局部 RT/ATAC。
- `sog1_targets`：旧分析中的 309 个 Arabidopsis SOG1 target，经 `rice2ara.csv` 转为 rice orthologues。
- `cell_cycle`：`Strictly_Matched_Cell_Cycle_Genes.csv` 中的 40 个水稻 cell-cycle genes。

## 输入

- 任务 009 gene-level RT/ATAC/RNA metrics：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/2026-08-26_009_cpp8-allele-path-replication/output/gene_level_path_metrics_009.tsv.gz`
- 任务 009 gene/TSS/RT mapping：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/2026-08-26_009_cpp8-allele-path-replication/output/gene_tss_rt_metrics_009.tsv.gz`
- DNA-repair panel、SOG1 targets、cell-cycle genes 和 rice–Arabidopsis mapping：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/`
- 参考基因组：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta`
- 注释：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3`

旧目录全部只读；本任务输出仅写入当前目录。

## 主要假设和终点

主 universe 为 TSS ±1 kb 内至少含一个 OCR、且通过任务 009 联合 RNA filter 的基因。MS 为预定义主阶段。

两个等位基因的保守共识扰动定义为：仅在 `oscpp8-1` 和 `oscpp8-3` 同号时取较小绝对效应，异号时记为 0。

主检验 family：

1. two-allele consensus TSS RT-instability magnitude；
2. two-allele consensus promoter MS accessibility-disruption magnitude；
3. 是否进入任务 012 定义的 top-5% RT-instability hotspot，阈值 `0.2632919444`。

共 3 个预定义基因集 × 3 个终点，使用 Benjamini–Hochberg FDR 校正。

## 统计模型

连续扰动 magnitude 经 `log1p` 和标准化后拟合：

```text
local disruption ~ gene-set membership
                 + WT WRT
                 + baseline RNA
                 + gene length
                 + promoter OCR number
                 + WT promoter ATAC
                 + chromosome
```

hotspot 采用相同协变量的 linear-probability model，系数解释为 adjusted risk difference。

- 主要区间和 P 值：2,000 次同步 1-Mb genomic-block bootstrap。
- 次级 matched null：在 chromosome × WT-WRT quintile × baseline-RNA quintile 内进行 5,000 次 label permutation。
- 随机种子：`20260826`。
- ±3 kb、ES/LS 和 allele-specific magnitude 为预定义敏感性分析。

## 运行

```bash
python code/01_damage_response_locality.py \
  --gene-metrics data/gene_level_path_metrics_009.tsv.gz \
  --gene-tss data/gene_tss_rt_metrics_009.tsv.gz \
  --dna-repair data/dna_repair_genes.csv \
  --cell-cycle data/cell_cycle_genes.csv \
  --sog1 data/sog1_targets_arabidopsis.csv \
  --rice2ara data/rice2ara.csv \
  --output-dir output \
  --plot-dir plots \
  --bootstrap 2000 \
  --permutations 5000
```

## 软件环境

- Python 3.11.4
- NumPy 1.26.4
- pandas 1.5.3
- SciPy 1.10.1
- statsmodels 0.14.0
- Matplotlib 3.7.1

## 主要结果

### 基因集和分析 universe

- 主 ±1 kb promoter universe：22,043 genes，379 个 1-Mb genomic blocks。
- ±3 kb sensitivity universe：23,460 genes。
- 主 universe 中的 gene-set size：
  - manuscript DNA-repair panel：9 genes；
  - rice orthologues of SOG1 targets：133 genes；
  - cell-cycle genes：39 genes。
- gene-set overlap 很小：DNA-repair panel 与 SOG1 targets 仅重叠 1 gene；SOG1 targets 与 cell-cycle genes 重叠 2 genes。

### 主检验

9 个主检验均未通过 FDR 校正，最小 `inference FDR=0.499`。

| gene set | endpoint | adjusted effect | 95% interval | inference P | FDR |
|---|---|---:|---:|---:|---:|
| DNA-repair panel | consensus RT instability | 0.334 s.d. | −0.254–0.810 | 0.222 | 0.499 |
| DNA-repair panel | consensus MS ATAC disruption | −0.078 s.d. | −0.840–0.719 | 0.794 | 0.939 |
| SOG1 targets | consensus RT instability | 0.104 s.d. | −0.040–0.248 | 0.157 | 0.499 |
| SOG1 targets | consensus MS ATAC disruption | −0.106 s.d. | −0.261–0.053 | 0.196 | 0.499 |
| cell-cycle genes | consensus RT instability | 0.167 s.d. | −0.124–0.461 | 0.283 | 0.509 |
| cell-cycle genes | consensus MS ATAC disruption | 0.123 s.d. | −0.159–0.426 | 0.398 | 0.597 |

DNA-repair panel 的 RT point estimate 为正，但仅有 9 genes，block-bootstrap interval 很宽且跨 0；不能作为局部富集证据。

### RT-instability hotspot

使用任务 012 的全基因组 top-5% threshold 后，主 promoter universe 中仅 85/22,043 genes（0.386%）落入 hotspot，三个 gene sets 均为 0 hit。由于这是 rare/zero-event endpoint，普通 block-linear model 会低估 label-assignment uncertainty。因此：

- 连续 RT/ATAC endpoint 以 1-Mb block bootstrap 为 primary inference；
- hotspot endpoint 以 covariate-matched label permutation 为 primary inference；
- block-linear hotspot interval 仅保留在完整表中，不用于显著性结论。

Matched P values：DNA-repair panel `P=0.214`、SOG1 targets `P=0.998`、cell-cycle genes `P=0.835`；均不显著。

进一步使用全基因组 top-10% threshold，或在 promoter universe 内重新定义 relative top-5%，仍无 gene set 通过 FDR（全部 `FDR≥0.409`）。

### 敏感性分析

- ±3 kb promoter window：无结果通过 FDR。
- ES、MS、LS two-allele consensus ATAC disruption：无阶段通过 FDR。
- 两个 CPP8 alleles 分别分析：无结果通过 family-wise BH correction。
- SOG1 targets 在 `oscpp8-3` 中有 nominal RT magnitude signal，但未通过 allele-sensitivity FDR，且缺少 `oscpp8-1` 重复，因此不作为结论。
- RNA magnitude 仅为 descriptive positive control。DNA-repair panel 是按 RNA outcome 选取的集合，不能把其 RNA enrichment 作为独立验证。

## 最终解释

本分析未发现 DNA-repair panel、SOG1 targets 或 cell-cycle genes 优先位于两个 CPP8 alleles 共同发生强 RT instability 或 promoter accessibility disruption 的位点。

因此，现有数据更支持以下有边界的解释：

> CPP8 loss produces reproducible local RT–ATAC coupling across the genome, but the DNA-damage and cell-cycle transcriptional programmes are not preferentially concentrated at the strongest local RT/ATAC perturbations. These programmes are therefore more consistent with a distributed or trans-acting response to replication stress than with uniform direct regulation through local RT shifts.

这里的 `more consistent with` 不是因果证明。bulk、非配对设计不能区分 checkpoint signalling、细胞状态变化和其他平行 CPP8 effects。

## 稿件位置

- 结果功能分类：`qualification`，会实质限制 Figure 8 的 direct/local mechanism claim。
- 正文：在 Figure 8 的 DNA-repair/cell-cycle Results 末尾保留一句限定。
- Discussion：用一句区分 genome-wide local RT–ATAC coupling 与 distributed stress-response programme。
- 完整 forest plot、stage heatmap、hotspot sensitivity 和全部统计表：Supplementary Information。
- 不建议占用新的主图 panel；该结果的价值是限定机制，而不是建立第二条正面主线。

英文 Results、Discussion、Methods 和图注候选见：
`output/manuscript_integration_014.md`。

## 解释边界

- 阳性结果表示预定义转录程序同时富集局部 RT/ATAC 扰动，支持局部耦联，但不能证明 RT 直接导致该基因的表达变化。
- 阴性结果表示 DNA-damage/cell-cycle RNA programme 并不集中于局部 RT/ATAC 强扰动基因，更符合全局或 trans-acting replication-stress response。
- bulk RT、ATAC 和 RNA 来自不同样本集合；1-Mb block inference 控制空间伪重复，但不能替代生物学重复。
- DNA-repair panel 只有 9 genes，统计区间较宽。
- SOG1 set 由 Arabidopsis targets 经本地 orthology table 映射到水稻，不能等同于已在水稻 root tip 中实验确认的 direct SOG1 targets。
- hotspot 在 promoter-gene universe 中稀少，结果主要依赖 matched permutation，而不是 zero-event linear-model interval。

## 输出

- `output/damage_response_locality_primary_014.tsv`
- `output/damage_response_locality_sensitivity_014.tsv`
- `output/damage_response_locality_models_014.tsv`
- `output/gene_universe_1000bp_014.tsv.gz`
- `output/gene_universe_3000bp_014.tsv.gz`
- `output/gene_set_membership_014.tsv`
- `output/gene_set_overlap_014.tsv`
- `output/damage_response_locality_summary_014.json`
- `plots/damage_response_locality_014.pdf|svg|png`
- `plots/damage_response_stage_heatmap_014.pdf|svg|png`
- `output/logs/02_damage_response_locality_final.log`
