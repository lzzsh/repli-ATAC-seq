# 008 — CPP8 RT–ATAC–RNA 主图与稿件整合

## 目的

将任务 003–006 的 CPP8 统计结果组织为一张承担单一机制主线的主文复合图，并明确该图在当前稿件 Results、Methods、Discussion、图注和补充材料中的插入位置与修改方式。

## 写作框架

本任务使用本地安装的 `nature-writing` 技能，按 Nature Portfolio 风格的“最短充分证据链”组织结果：

```text
可重复遗传表型
→ 同位点阶段化 RT–ATAC 耦联
→ promoter RT–ATAC–RNA 统计路径
→ 跨染色体预测与空间负对照
→ 明确因果边界
```

## 输入

- 稿件：`data/manuscript_20260705.pdf`
- 任务 003 图：`data/source_003_mutant_rt_robustness.pdf`
- 任务 004 图：`data/source_004_rt_atac_stage_association.pdf`
- 任务 005 图：`data/source_005_rt_atac_rna_paths.pdf`
- 任务 006 图：`data/source_006_causal_path_validation.pdf`
- 详细统计报告：项目根目录 `CPP8_RT_ATAC_RNA_DETAILED_REPORT.md`

以上 `data/` 文件均为指向既有只读结果的软链接；本任务未修改来源文件。

## 主要输出

- `MANUSCRIPT_FIGURE_INTEGRATION_PLAN.md`：总图设计、正文插入位置、英文 Results、图注、Methods、Discussion和补充材料分配。
- `plots/previews/`：四张来源图的PNG预览，仅用于设计检查。

## 核心决定

- 推荐将总图设为重新设计后的 **Figure 7**，而不是新增第九张主图。
- 当前 Figure 7 中 DAR 数量、Venn、注释、普通散点等描述性面板移至补充材料。
- 新 Figure 7 使用 7 个面板，聚焦 CPP8；CPP11异质性不与CPP8合并解释。
- 主文报告决定性效应与主要区间；完整bootstrap、置换、窗口敏感性和模型诊断放补充材料或Source Data。
- 采用“genetic perturbation–anchored statistical path”表述，不称为严格 causal mediation。

## 未执行内容

本任务仅完成图文整合规划，尚未重新绘制或排版最终复合图。后续若制作正式Figure 7，应在本目录的 `code/`、`plots/` 和 `output/logs/` 中完成，不能覆盖任务003–006的来源图。

## 权限

目录及所有实体文件保持 owner-only；group/other 无读、写、执行权限，子目录设置私有默认ACL。
