# 013 — RT–ATAC–RNA路径方向性的反证与兼容性检验

## 状态

**正式分析已完成（2026-08-26）。** 本任务检验候选排序的非对称兼容性，不从非配对bulk组学中识别真实因果DAG。

## 科学问题与结论

比较三种路径：

```text
F:  ΔRT → promoter ΔATAC → RNA
R1: promoter ΔATAC → ΔRT → RNA
R2: RNA → promoter ΔATAC → ΔRT
```

MS中forward路径在两个alleles和两个promoter窗口均稳定，并显著强于ATAC-first排序；但RNA-first路径本身也稳定且接近forward。结合LOCO、时间/空间结构和负对照，数据更兼容具有局部、MS阶段结构的`RT–ATAC–RNA` ordering，但不能排除RNA-first统计排序或CPP8缺失对三层的平行影响。

## 输入

- gene metrics：`../2026-08-26_009_cpp8-allele-path-replication/output/gene_level_path_metrics_009.tsv.gz`
- stage contrasts：`../2026-08-26_010_rt-stage-match-coupling/output/stage_match_primary_contrasts_010.tsv`
- scale contrasts：`../2026-08-26_011_rt-atac-spatial-scale/output/multiscale_primary_contrasts_011.tsv`
- offset curves：`../2026-08-26_011_rt-atac-spatial-scale/output/spatial_offset_curve_011.tsv.gz`

## 模型与推断

每个ordering分别拟合两段标准化线性路径，控制WT WRT、WT ATAC、WT RNA、启动子OCR数量和染色体。三种ordering在同一次2,000次1-Mb block bootstrap中估计；379个blocks；seed=`20260826`。比较`IF−IR1`和`IF−IR2`的同步bootstrap CI，而不是比较显著与不显著。

LOCO按12条染色体留一预测RNA或RT outcome。负对照使用wrong-stage mediator以及染色体内circular-shift RNA；每项最多1,000次bootstrap。

## 正式结果

### MS ±1 kb路径

| allele | ordering | product path | 95% block-bootstrap CI |
|---|---|---:|---:|
| `oscpp8-1` | RT→ATAC→RNA | 0.01628 | 0.01276–0.02017 |
| `oscpp8-1` | ATAC→RT→RNA | −0.00381 | −0.00530至−0.00243 |
| `oscpp8-1` | RNA→ATAC→RT | 0.01481 | 0.01168–0.01806 |
| `oscpp8-3` | RT→ATAC→RNA | 0.02656 | 0.02226–0.03128 |
| `oscpp8-3` | ATAC→RT→RNA | −0.00424 | −0.00616至−0.00224 |
| `oscpp8-3` | RNA→ATAC→RT | 0.02280 | 0.01907–0.02638 |

±3 kb保持相同格局：forward分别为0.01833和0.02766；RNA-first分别为0.01643和0.02282。

### Forward−reverse contrasts

- MS ±1 kb，forward−ATAC-first：
  - `oscpp8-1` 0.02010，95% CI 0.01574–0.02481；
  - `oscpp8-3` 0.03080，95% CI 0.02583–0.03642。
- MS ±1 kb，forward−RNA-first：
  - `oscpp8-1` 0.00147，95% CI 0.00023–0.00294；
  - `oscpp8-3` 0.00377，95% CI 0.00200–0.00590。

forward虽略强于RNA-first，但差距远小于相对ATAC-first的差距；不能据此证明时间方向。

### 条件LOCO预测

MS ±1 kb RNA outcome：

- 已知ΔRT后加入ΔATAC：mean MSE reduction为`oscpp8-1=0.01045`、`oscpp8-3=0.01756`；paired Wilcoxon分别`P=0.000244`和`0.000488`。
- 已知ΔATAC后加入ΔRT：0.000744和0.000478；两者`P=0.0461`，但增量约小一个数量级。
- pooled LOCO中，ΔATAC相对baseline的ΔR²为0.02165和0.03653；ΔRT相对baseline为0.00056和−0.00007。

这支持ATAC是RNA预测的更直接层，但不构成方向识别。

### 时间、空间与负对照

- 两alleles均有global MS dominance；没有共同diagonal matching。
- 六个offset曲线均在0 kb达峰，支持locus alignment。
- circular-shift RNA后MS路径均崩溃且CI跨0。
- wrong-stage mediator明显减弱forward path，但`oscpp8-3`仍有残余正路径，与任务010的全局MS/跨阶段耦联一致，而非绝对stage specificity。

## 运行命令

```bash
python3 code/02_fit_directional_models.py \
  --gene-metrics ../2026-08-26_009_cpp8-allele-path-replication/output/gene_level_path_metrics_009.tsv.gz \
  --stage-contrasts ../2026-08-26_010_rt-stage-match-coupling/output/stage_match_primary_contrasts_010.tsv \
  --scale-contrasts ../2026-08-26_011_rt-atac-spatial-scale/output/multiscale_primary_contrasts_011.tsv \
  --offsets ../2026-08-26_011_rt-atac-spatial-scale/output/spatial_offset_curve_011.tsv.gz \
  --output-dir output --plot-dir plots --bootstrap 2000
```

## 软件环境

Python 3.11.4；NumPy 1.26.4；pandas 1.5.3；SciPy 1.10.1；scikit-learn 1.3.0；matplotlib 3.7.1。

## 稿件位置

- 完整forward/reverse forest、LOCO和negative controls放Supplementary Figure/Methods。
- 正文最多1–2句限制性总结：observed ordering比ATAC-first更兼容，但RNA-first和parallel perturbation不能排除。
- 这一结论改变主文措辞边界，因此Discussion必须保留，不可完全藏入SI。

## 禁止表述

- `causal direction proven`；
- `reverse causation excluded`；
- `CPP8 genotype is a valid instrument`；
- 用SEM/AIC或不同outcome的R²直接判定真实时间方向。

## 输出

- 表：`all_direction_path_effects_013.tsv`、`forward_vs_reverse_contrasts_013.tsv`、`directional_loco_summary_013.tsv`、`directional_loco_comparisons_013.tsv`、`negative_control_paths_013.tsv`、`directional_evidence_matrix_013.tsv`。
- 图：`forward_reverse_path_forest_013.pdf`。
- `02_fit_directional_models_smoke.log`为调试记录；正式数值以无`_smoke`日志和当前TSV为准。
