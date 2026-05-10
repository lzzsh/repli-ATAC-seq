# Basenji2 Repli-seq

基于 DNA 序列的 Basenji2 卷积模型，预测 Repli-seq ES/MS/LS/G1 信号（log1p TPM），并由预测值推导 WRT（Weighted Replication Timing）。架构对齐 Basenji2：one-hot 输入 → conv stem → conv tower → dilated residual tower → 中心 8 bin 逐位置线性预测头。

每个训练样本以一个 GFF3 标注 bin 为中心截取 32 kb 窗口，模型输出中心 8 个 128 bp bin（共 1024 bp）的预测值，loss 在这 8 个 bin 上计算。

---

## 目录

1. [环境安装](#1-环境安装)
2. [数据准备](#2-数据准备)
3. [配置文件](#3-配置文件)
4. [训练](#4-训练)
5. [评估](#5-评估)
6. [项目结构](#6-项目结构)
7. [模型框架](#7-模型框架)

---

## 1. 环境安装

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install pyfaidx pandas numpy scipy scikit-learn pyyaml tqdm tensorboard
```

Python 3.10+，PyTorch 2.x。

---

## 2. 数据准备

### 2.1 参考基因组

将各物种的参考基因组 FASTA 文件放到对应目录，并建立 pyfaidx 索引：

```bash
python -c "from pyfaidx import Fasta; Fasta('/path/to/genome.fa')"
```

### 2.2 Repli-seq count 文件

每个物种需要一个 TSV 文件，每行对应一个 128 bp bin，格式如下：

```
#chrom  start   end     ES_count    MS_count    LS_count    G1_count
chr01   0       128     148         150         368         150
chr01   128     256     680         1945        2987        1145
...
```

- `chrom`：染色体名称，需与 FASTA 中的 sequence name 一致（列名 `#chrom` 会自动重命名）
- `start` / `end`：bin 的基因组坐标（0-based，半开区间），bin size = 128 bp
- `ES_count` / `MS_count` / `LS_count` / `G1_count`：各 phase 的 raw read count

> **归一化由代码自动完成**：TPM 归一化 → log1p 变换 → WRT 计算。GFF3 中标注为边界（ESMS/MSLS/ESLS/ESMSLS）的 bins 会被自动过滤。

### 2.3 数据划分说明

划分以**完整染色体**为单位，禁止随机划分（相邻 bins 序列高度相关，随机划分会造成 genomic leakage）。

默认划分（见 `data/manifest.yaml`）：

| 物种 | train | val | test |
|------|-------|-----|------|
| 水稻 | chr01–chr08 | chr09–chr10 | chr11–chr12 |

---

## 3. 配置文件

### 3.1 manifest（物种清单）

`data/manifest.yaml` 定义每个物种的数据路径和染色体划分：

```yaml
species:
  - name: rice
    fasta: /path/to/genome.fasta
    count_tsv: data/labels/rice_wt_counts.tsv
    gff3: /path/to/replication_phase.gff3
    train_chroms: [chr01, chr02, chr03, chr04, chr05, chr06, chr07, chr08]
    val_chroms: [chr09, chr10]
    test_chroms: [chr11, chr12]
```

### 3.2 训练配置

`src/configs/transformer_wt.yaml` 包含所有超参数，关键字段：

```yaml
data:
  input_window_length: 32768   # 32 kb 输入窗口，对齐 Basenji2

augmentation:
  rc_prob: 0.5                 # reverse complement 概率
  shift_max: 3                 # 训练时随机平移最大 bins 数

model:
  bn_momentum: 0.1             # TF bn_momentum=0.9 → PyTorch 1-0.9=0.1

loss:
  lambda_phase: 1.0            # MSE loss 权重

training:
  learning_rate: 0.15          # SGD lr，对齐 Basenji2
  momentum: 0.99
  batch_size: 8
  early_stopping_patience: 16  # val_loss 连续 16 个 epoch 不下降则停止
  mixed_precision: true
```

---

## 4. 训练

### 4.1 单卡训练

```bash
python scripts/train.py --config src/configs/transformer_wt.yaml
```

### 4.2 多卡 DDP 训练

```bash
torchrun --nproc_per_node=4 scripts/train.py --config src/configs/transformer_wt.yaml
```

### 4.3 断点续训

```bash
torchrun --nproc_per_node=4 scripts/train.py \
    --config src/configs/transformer_wt.yaml \
    --resume outputs/basenji2_wt/checkpoints/best_model.pt
```

### 4.4 训练过程监控

```bash
tensorboard --logdir logs/
```

TensorBoard 记录：
- `train/loss_total`、`train/loss_phase`（epoch 平均）
- `train/lr`
- `val/loss_total`、`val/loss_phase`
- `val/phase_pearson_ES`、`val/phase_pearson_MS`、`val/phase_pearson_LS`、`val/phase_pearson_G1`
- `val/mean_phase_pearson`、`val/wrt_pearson`

### 4.5 checkpoint

最优模型保存在 `outputs/basenji2_wt/checkpoints/best_model.pt`，由 `val_loss_total`（越低越好）决定。

---

## 5. 评估

### 5.1 WT 测试集评估

```bash
python scripts/evaluate.py \
    --mode wt \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/wt_test_metrics.tsv`，包含：
- `phase_pearson_ES/MS/LS/G1`、`phase_spearman_ES/MS/LS/G1`
- `mean_phase_pearson`、`mean_phase_spearman`
- `wrt_pearson`、`wrt_spearman`

### 5.2 Leave-one-species-out 评估

```bash
python scripts/evaluate.py \
    --mode cross_species \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --held_out rice \
    --output outputs/metrics
```

---

## 6. 项目结构

```
src/
  data/
    data_utils.py       # GenomeSequence, one_hot_encode, TPM norm, WRT, labels
    dataset.py          # SpeciesConfig/Manifest, RepliSeqDataset
  models/
    model.py            # Basenji2Model, _Basenji2Trunk, _DilatedResidual, _SharedHead, PhaseLoss
  configs/
    transformer_wt.yaml # 训练超参数
  trainer.py            # DDP 训练循环 + TensorBoard
  eval.py               # 评估指标 + WT/cross-species 评估

scripts/
  train.py              # 训练入口（支持 torchrun）
  evaluate.py           # 评估入口

data/
  manifest.yaml         # 物种清单（路径 + 染色体划分）
  genomes/              # 参考基因组 FASTA（不入 git）
  labels/               # Repli-seq count TSV（不入 git）

outputs/                # 模型 checkpoint、评估结果（不入 git）
logs/                   # TensorBoard 日志（不入 git）
```

---

## 7. 模型框架

### 7.1 整体结构

```
DNA sequence [32 kb]，以 GFF3 标注 bin 为中心截取
      │
 one_hot_encode
 [B, 4, 32768]  (A/C/G/T，N→全零)
      │
 Conv stem
 Conv1d(4→288, k=15, pool=2)          → [B, 288, 16384]
      │
 Conv tower（6 层，每层 pool=2）
 filters: 288→339→399→470→554→652→768  (×1.1776 递增)
 每层: BN → GELU → Conv1d(k=5, pool=2)
                                       → [B, 768, 256]  @ 128bp/bin
      │
 Dilated residual tower（11 层）
 bottleneck 768→384→768，dilation rate_mult=1.5
                                       → [B, 768, 256]
      │
 Cropping1D（两端各裁 16 bins）        → [B, 768, 224]  @ 128bp/bin
      │
 Bottleneck conv: BN→GELU→Conv1d(768→1536, k=1, dropout=0.05)
                                       → [B, 1536, 224]
      │
 取中心 8 个 bin（bin 108:116）        → [B, 1536, 8]
 对应基因组位置: win_center ± 512bp，共 1024bp
      │
 Linear(1536→4) 逐 bin 独立应用
      │
 phase_pred [B, 8, 4]   (ES/MS/LS/G1 log1p TPM)
      │
 loss = MSE(phase_pred, phase_labels[B, 8, 4])
```

### 7.2 采样策略

每个训练样本以一个 GFF3 标注 bin 为中心截取 32 kb 窗口：

```
标注 bin（ES/MS/LS）:  ████  ← 窗口中心 = 模型感受野中心
窗口:  |←────────16384bp────────→|████|←────────16384bp────────→|
标签:  中心 8 个 128bp bin 的 [ES, MS, LS, G1] log1p
```

- 样本数 = GFF3 标注 bin 数（边界 bin 已过滤）
- 标注 bin 严格落在中心 8 个输出 bin 内，无坐标偏移

### 7.3 Dilated Residual Block

```
x [B, 768, L]
│
├─ GELU → Conv(k=3, dilation=d) → BN   ← 扩大感受野
│         │
│  GELU → Conv(k=1, dropout=0.3) → BN  ← gamma 零初始化（InitZero）
│
x = x + residual
```

11 层 dilation 序列（rate_mult=1.5，取整）：1, 2, 2, 3, 5, 8, 11, 17, 25, 38, 57，最大感受野约 ±22 kb。

### 7.4 PhaseLoss

```
loss = MSE(phase_pred, phase_labels)
```

`phase_labels` shape `[B, 8, 4]`，为 log1p(TPM) 归一化后的 ES/MS/LS/G1 值，8 个 bin 等权。

### 7.5 训练策略（对齐 Basenji2）

| 参数 | 值 |
|------|----|
| Optimizer | SGD + momentum |
| lr | 0.15 |
| momentum | 0.99 |
| Scheduler | ReduceLROnPlateau（patience=16，factor=0.2） |
| Clip norm | 2.0 |
| Batch size | 8 |
| Early stopping patience | 16 |
| 数据增强 | reverse complement（p=0.5）+ random shift（±3 bins） |
