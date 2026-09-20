# 010 — WT复制时期与ATAC测量阶段的匹配及滞后分析

## 状态

**正式分析已完成（2026-08-26）。** 本任务只检验时间结构，不重复总体RT极化或promoter路径分析。

## 科学问题

用`WT RT class × ATAC stage`的3×3设计区分三种解释：局部复制阶段同步、复制后一阶段滞后、或不依赖局部RT类别的全局MS敏感性。

## 输入与数据结构

- 主输入：`../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz`
- 每个OCR按WT ES/MS/LS S/G1 normalized signal的最大值分为E、M或L；confidence为最大阶段比例减次大阶段比例。
- 可用OCR的WT类别计数约为E 12.9k、M 35.7k、L 40.2k；每个等位基因按三个ATAC阶段堆叠。

## 模型与推断

```text
z(ΔATAC) ~ 0 + RTclass×ATACstage cell slopes for z(ΔRT)
           + WT WRT + stage-specific WT ATAC
           + RT confidence + log(OCR length+1)
           + chromosome fixed effects
```

- 3×3 cell-specific ΔRT slopes分别在两个alleles中估计。
- 预定义contrast：diagonal dominance、+1-stage lag、global MS dominance。
- 2,000次同步1-Mb block bootstrap；379个blocks；seed=`20260826`。
- 六个allele-specific primary contrasts进行BH FDR；三个allele-difference contrasts单独校正。

## 正式结果

| allele | contrast | estimate | 95% block-bootstrap CI | FDR |
|---|---|---:|---:|---:|
| `oscpp8-1` | diagonal dominance | −0.0041 | −0.0150–0.0068 | 0.573 |
| `oscpp8-3` | diagonal dominance | −0.0376 | −0.0495至−0.0259 | 0.00150 |
| `oscpp8-1` | +1-stage lag | 0.0016 | −0.0137–0.0170 | 0.875 |
| `oscpp8-3` | +1-stage lag | 0.0449 | 0.0283–0.0615 | 0.00150 |
| `oscpp8-1` | global MS dominance | 0.0219 | 0.0120–0.0324 | 0.00150 |
| `oscpp8-3` | global MS dominance | 0.0895 | 0.0769–0.1027 | 0.00150 |

两个alleles的global MS dominance差值为−0.0677（`oscpp8-1−oscpp8-3`），95% CI −0.0833至−0.0523，说明共同方向成立但强度allele-dependent。

3×3 cell结果进一步显示：

- `oscpp8-3`的最大cell slope位于L-class OCR在MS测量时，β=0.2110（95% CI 0.1919–0.2303），而非对角线L–LS。
- `oscpp8-1`中L–MS也较强，β=0.1166（0.0995–0.1335）。
- M-class在`oscpp8-1`三个阶段的slopes均接近0，削弱统一local-stage matching解释。

## 结论

两个等位基因共同支持“MS是CPP8缺失后RT–ATAC耦联的全局敏感阶段”。数据不支持跨allele一致的local RT-class diagonal matching；+1-stage lag只在`oscpp8-3`成立，不能写成CPP8通用的复制后染色质恢复机制。

## 运行顺序

```bash
python3 code/01_build_stage_match_matrix.py \
  --ocr-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz \
  --output output/stage_match_matrix_010.tsv.gz
python3 code/02_fit_stage_match_models.py \
  --matrix output/stage_match_matrix_010.tsv.gz \
  --output-dir output --bootstrap 2000
python3 code/03_plot_stage_match_results.py \
  --matrix output/stage_match_matrix_010.tsv.gz \
  --coefficients output/stage_match_coefficients_010.tsv \
  --contrasts output/stage_match_primary_contrasts_010.tsv \
  --lag output/stage_lag_distance_010.tsv \
  --output-dir output --plot-dir plots
```

## 软件环境

Python 3.11.4；NumPy 1.26.4；pandas 1.5.3；SciPy 1.10.1；statsmodels 0.14.0；matplotlib 3.7.1；seaborn 0.12.2。

## 稿件位置与边界

- 正文或主Figure 7只保留global MS dominance及其双allele一致性。
- 完整3×3 heatmap、lag-distance和trajectory进入Supplementary Figure/Source Data。
- 独立阶段样本不是同一细胞longitudinal trajectory；“MS最强”也不能写成已证明局部复制后的时间因果顺序。

## 输出

- 表：`stage_match_coefficients_010.tsv`、`stage_match_primary_contrasts_010.tsv`、`stage_lag_distance_010.tsv`、`stage_match_allele_heterogeneity_010.tsv`、`delta_atac_trajectories_010.tsv`。
- 图：`rt_class_by_atac_stage_heatmap_010.pdf`、`stage_lag_contrasts_010.pdf`、`delta_atac_trajectories_010.pdf`。
- `02_fit_stage_match_models_smoke.log`为调试记录；正式数值以无`_smoke`日志和当前TSV为准。
