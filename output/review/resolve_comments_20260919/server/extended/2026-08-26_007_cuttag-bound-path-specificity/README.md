# 007 — CUT&Tag-bound promoter 路径特异性

## 目的

检验 CPP8/CPP11 CUT&Tag-bound promoters 的 RT→ATAC→RNA 间接路径是否显著强于未结合 promoters。

## 输入

- 任务 005 gene-level metrics。
- CPP8 peaks：`cuttag_tcx/macs2/macs2_p1e-5/tcx_cut_out/tcx_cut_peaks_sort.narrowPeak`
- CPP11 peaks：`cuttag_tcx/macs2/macs2_p1e-5/sol_cut_out/sol_cut_peaks_sort.narrowPeak`

## 方法

- CUT&Tag peak 与 TSS ±1/±3 kb 任一重叠定义为 bound promoter。
- a 和 b 路径中加入 `bound × ΔRT`、`bound × ΔATAC` interaction。
- 同步 1-Mb block bootstrap 估计 bound 与 unbound 的间接效应及其差值。

## 主要结果

- ±1 kb 中，CPP8 和 CPP11 分别有约 1,118 和 1,236 个 bound expressed promoters；±3 kb 分别约 2,135 和 2,309 个。
- 12 个组合中没有一个 `bound − unbound` 间接效应在多重校正后显著；所有 contrast FDR ≥ 0.448。
- CPP8 ES/MS 的 bound path 数值上略大，但 block-bootstrap CI 跨 0。

## 解释与限制

当前结果不支持把全体 CPP CUT&Tag-bound promoter 写成统一增强的 causal-path subset。可能原因包括：

- CUT&Tag 来自 protoplast，而 RT/ATAC/RNA 来自 root-tip，细胞背景不匹配。
- binding 本身不是充分条件，方向和效应可能依赖位点上下文。
- bound/unbound 交互检验比只在 bound subset 内报告相关性更严格。

本结果适合作为负对照或限制，不建议作为主图正面结论。若未来获得 root-tip CUT&Tag/CUT&RUN，再重新检验直接靶位点特异性。

## 输出

- `output/cuttag_bound_promoter_map_007.tsv.gz`
- `output/cuttag_bound_path_specificity_007.tsv`
- `plots/cuttag_bound_path_specificity_007.pdf`

