# 016 — Paper narrative logic audit

## 目的

通读 `Repli-ATAC-seq-20260705.pdf` 全文（含 Abstract、Introduction、Results 的 Figure 1–8 各小节及每张大图和小图、Discussion、Methods、Figure S1–S5），梳理全文行文逻辑，逐图核对其在证据链中的位置，并评估逻辑是否正确。

本任务只做文本与已有审计文档的整合判断，不调用视觉模型。

## 输入

- 稿件：`/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/Repli-ATAC-seq-20260705.pdf`
- 正文文本：`output/paper_text.txt`（pdftotext -layout 提取，1244 行）
- 已有审计：`2026-08-27_015_manuscript-figure-audit/output/MANUSCRIPT_FIGURE_AUDIT_015.md`
- 已有因果路径报告：`CPP8_RT_ATAC_RNA_DETAILED_REPORT.md`、`MANUSCRIPT_STRENGTHENING_SUMMARY.md`

## 输出

- `output/paper_logic_assessment_016.md`：全文行文逻辑梳理 + 逐图证据链定位 + 逻辑正确性评估 + 逻辑断裂点清单。
- `output/cross_reference_check_016.tsv`：正文图引用与图注一致性抽检。
- `output/point4_point6_detail_016.md`：第 4、6 点详细展开。
- `output_v2/s4s5_fix_verification_016.md`：新版 PDF 中 S4/S5 引用修正的 OCR 核验记录。
- `output_v2/point4_revision_016.md`：第 4 点采纳作者对 8F 图注限定的更正后的修订。

## 新版 PDF 更新（2026-08-27）

作者重导 PDF：49 页、全图栅格化（无文本层，需 OCR）。已核验：第 6 点（S4/S5）已修正并通过 OCR 确认；第 4 点已采纳作者关于 8F 图注范围正确限定（OsCPP8-associated promoter-proximal）的更正，撤回“global 全基因组”表述，其余证据强度论据保留。
