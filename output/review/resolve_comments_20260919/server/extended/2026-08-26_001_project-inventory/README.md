# 001 Repli-ATAC-seq project inventory

## 任务目的

解析 `extended` 中的 Repli-ATAC-seq 论文 PDF，盘点同级 `Projects` 路径下与论文相关的既有数据、代码和结果，并建立可继续维护的项目级 Markdown 档案。

## 创建日期

2026-08-26（Asia/Shanghai）

## 输入

- 原始论文 PDF：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/Repli-ATAC-seq-20260705.pdf`
- 本任务软链接：
  `data/Repli-ATAC-seq-20260705.pdf`
- 既有分析来源根目录（本次仅只读检查）：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects`
- 目录和保存规范参考：
  `/storage2/liuxiaodongLab/liaozizhuo/Projects/BowenHu/hypoblast_rebuttal/AGENTS.md`

## 本次检查方式

- 用 `pdfinfo` 检查 PDF 元数据、页数和文件状态。
- 用 `pdftotext` 提取题名、摘要、图注、方法、软件版本和数据线索。
- 用 `find`、`du`、`rg`、`sed` 检查目录体量、文件数量、脚本、样本表、结果文件和跨目录依赖。
- 用 `sha256sum` 记录 PDF 校验值。
- 未运行任何旧分析脚本，未修改任何既有项目目录。

## 交付物

- 项目工作规范：`../AGENTS.md`
- 项目数据与分析索引：`../DATA_INVENTORY.md`
- 原始 PDF 的只读入口：`data/Repli-ATAC-seq-20260705.pdf`

## 当前目录结构

```text
2026-08-26_001_project-inventory/
├── README.md
├── code/
├── data/
│   └── Repli-ATAC-seq-20260705.pdf -> ../../Repli-ATAC-seq-20260705.pdf
├── output/
│   └── logs/
└── plots/
```

## 主要结论

- 已定位论文主数据链：Repli-ATAC-seq、突变体 Repli-seq、阶段化 ATAC-seq、RNA-seq、CPP8/CPP11 与 WOX11 CUT&Tag，以及少量 Repliformer motif/ISM 后处理材料。
- 已从代码确认旧命名 `TCX2 = OsCPP8`、`SOL1 = OsCPP11`，并整理主要样本映射。
- 已发现若干复现风险，包括旧脚本硬编码输出路径、部分零字节结果、论文与当前代码的 DEG 阈值不一致，以及 `repliformer` 目录缺少完整训练代码和模型权重。
- 详细结果见项目根目录 `DATA_INVENTORY.md`。

## 后续更新规则

后续若补充样本表、图源、模型代码或 accession，优先更新 `DATA_INVENTORY.md`，并在新的规范任务目录中保存验证脚本和派生结果。

## 权限更新

2026-08-26 将本次确认的全部 Repli-ATAC-seq 相关目录递归设为所有者私有：移除 group/other 的全部权限，清除旧扩展 ACL，并为所有现有子目录设置私有默认 ACL。随后补查并纳入此前遗漏的 `GWAS` 目录。
