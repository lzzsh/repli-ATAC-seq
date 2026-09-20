# 011 — CPP8 RT–ATAC耦联的空间尺度与距离衰减

## 状态

**正式分析已完成（2026-08-26）。** 本任务只回答局部与邻域尺度，不重复时间匹配或基因组特征富集。

## 科学问题与结论

检验CPP8相关ΔATAC更接近同位点ΔRT，还是5–100 kb邻域或候选RT边界尺度的ΔRT。

正式结果支持局部耦联：六个allele×stage组合的offset曲线均在0 kb达到峰值；1-kb效应通常高于25/100-kb尺度；在同时控制四层邻域ΔRT后，local coefficient六项均为正且CI排除0。候选RT boundary附近未出现跨allele稳定增强。

## 输入

- OCR metrics：`../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz`
- RT bins：`../2026-08-26_009_cpp8-allele-path-replication/output/repli_rt_bin_metrics_009.tsv.gz`

## 空间变量

- Nested radii：1、5、25、100 kb。
- Non-overlapping components：local、1–5、5–20、20–50、50–100 kb。
- Offset curve：−100至+100 kb，5-kb步长；每点使用5-kb窗口。
- Candidate RT boundary proxy：WT WRT的20-kb rolling局部一阶差分极大值；primary阈值top 10%，敏感性top 5%和15%。near定义≤10 kb，interior定义≥50 kb。

boundary仅是连续RT信号的计算代理，不是Hi-C/TAD/compartment boundary。

## 模型与推断

```text
z(ΔATACstage) ~ z(ΔRTwindow) + WT WRTwindow
                + WT stage ATAC + log(OCR length+1)
                + chromosome fixed effects
```

同时模型将五个非重叠空间成分一起放入。所有主要尺度在同一次1-Mb block bootstrap中同步重采样，以保留尺度间协方差。正式bootstrap 1,000次；379个blocks；seed=`20260826`。

## 正式结果

### 1-kb相对25/100-kb的local dominance

| allele | stage | contrast | 95% block-bootstrap CI |
|---|---|---:|---:|
| `oscpp8-1` | ES | 0.0478 | 0.0277–0.0682 |
| `oscpp8-1` | MS | 0.0297 | 0.0097–0.0495 |
| `oscpp8-1` | LS | 0.0664 | 0.0435–0.0894 |
| `oscpp8-3` | ES | 0.0508 | 0.0300–0.0714 |
| `oscpp8-3` | MS | 0.0228 | −0.0025–0.0451 |
| `oscpp8-3` | LS | 0.0825 | 0.0574–0.1060 |

除`oscpp8-3` MS外，五项CI排除0；该例仍同方向。每log10半径的β趋势在另外五项均显著为负，`oscpp8-3` MS为−0.0145（CI −0.0296–0.0020）。

### 同时local-vs-annulus模型

| allele | ES local β | MS local β | LS local β |
|---|---:|---:|---:|
| `oscpp8-1` | 0.0885 (0.0701–0.1071) | 0.1104 (0.0944–0.1257) | 0.1696 (0.1507–0.1869) |
| `oscpp8-3` | 0.0839 (0.0637–0.1025) | 0.1216 (0.1060–0.1396) | 0.1684 (0.1491–0.1897) |

六项local系数均通过bootstrap。部分邻域annulus系数有正有负，提示邻域结构复杂，但不改变local component稳定存在的结论。

### Offset与boundary

- 六条offset曲线峰值全部在0 kb。
- top-10% boundary proxy中，六个`near−interior` contrasts的95% CI均跨0。
- top-5%/10%/15%阈值未形成跨allele、跨stage一致的boundary增强模式。

## 运行顺序

```bash
python3 code/01_build_multiscale_rt_metrics.py \
  --ocr-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/ocr_rt_atac_metrics_009.tsv.gz \
  --rt-bins ../2026-08-26_009_cpp8-allele-path-replication/output/repli_rt_bin_metrics_009.tsv.gz \
  --output output/multiscale_rt_metrics_011.tsv.gz
python3 code/02_fit_spatial_scale_models.py \
  --metrics output/multiscale_rt_metrics_011.tsv.gz \
  --output-dir output --bootstrap 1000
python3 code/03_rt_boundary_analysis.py \
  --metrics output/multiscale_rt_metrics_011.tsv.gz \
  --rt-bins ../2026-08-26_009_cpp8-allele-path-replication/output/repli_rt_bin_metrics_009.tsv.gz \
  --output-dir output --bootstrap 1000
python3 code/04_plot_spatial_scale.py \
  --coefficients output/multiscale_rt_atac_coefficients_011.tsv \
  --offsets output/spatial_offset_curve_011.tsv.gz \
  --boundaries output/rt_boundary_sensitivity_011.tsv \
  --plot-dir plots
```

## 软件环境

Python 3.11.4；NumPy 1.26.4；pandas 1.5.3；SciPy 1.10.1；statsmodels 0.14.0；matplotlib 3.7.1；seaborn 0.12.2。

## 稿件位置与边界

- 主Figure 7可加入一个scale-decay或0-kb peak panel，支持“spatially concentrated at local regulatory regions”。
- annulus全模型、41点offset曲线和boundary阈值敏感性进入SI。
- 不能写成局部RT改变在单细胞中传播到ATAC，也不能把boundary proxy称作TAD边界。

## 输出

- 表：`multiscale_rt_atac_coefficients_011.tsv`、`multiscale_primary_contrasts_011.tsv`、`local_vs_annulus_models_011.tsv`、`spatial_offset_curve_011.tsv.gz`、`rt_boundary_sensitivity_011.tsv`、`candidate_rt_boundaries_011.tsv.gz`。
- 图：`rt_atac_scale_decay_011.pdf`、`rt_atac_offset_curve_011.pdf`、`rt_boundary_coupling_011.pdf`。
- `*_smoke.log`为调试记录；正式数值以无`_smoke`日志和当前TSV为准。
