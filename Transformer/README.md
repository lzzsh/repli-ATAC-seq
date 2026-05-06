# Basenji2 Repli-seq

基于 DNA 序列的 Basenji2 卷积模型，预测 Repli-seq ES/MS/LS 信号和 RT class，支持 in silico mutagenesis 解释序列特征对复制时序的贡献。架构对齐 Basenji2：one-hot 输入 → conv stem → conv tower → dilated residual tower → 多任务预测头。

---

## 目录

1. [环境安装](#1-环境安装)
2. [数据准备](#2-数据准备)
3. [配置文件](#3-配置文件)
4. [训练](#4-训练)
5. [评估](#5-评估)
6. [In silico mutagenesis](#6-in-silico-mutagenesis)
7. [基线模型](#7-基线模型)
8. [项目结构](#8-项目结构)
9. [模型框架](#9-模型框架)

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

每个物种需要一个 TSV 文件，每行对应一个 bin，格式如下：

```
chrom   start   end     ES_count    MS_count    LS_count
chr1    0       1000    142         87          23
chr1    1000    2000    98          201         56
...
```

- `chrom`：染色体名称，需与 FASTA 中的 sequence name 一致
- `start` / `end`：bin 的基因组坐标（0-based，半开区间）
- `ES_count` / `MS_count` / `LS_count`：各 phase 的 raw read count

> **归一化由代码自动完成**：TPM 归一化 → log1p 变换 → WRT 计算 → RT class 分配（ES/MS/LS）。EM/ML 边界 bins 会被自动过滤。

### 2.3 数据划分说明

划分以**完整染色体**为单位，禁止随机划分（相邻 bins 序列高度相关，随机划分会造成 genomic leakage）。

默认划分（见 `data/manifest.yaml`）：

| 物种 | train | val | test |
|------|-------|-----|------|
| 水稻 | chr1–chr8 | chr9–chr10 | chr11–chr12 |

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

model:
  bn_momentum: 0.1
  head_hidden_dim: 128
  head_dropout: 0.1

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

### 4.3 训练过程监控

```bash
tensorboard --logdir logs/
```

TensorBoard 记录：
- `train/loss_total`、`train/loss_phase`、`train/loss_class`（epoch 平均）
- `train/lr`
- `val/loss_total`、`val/loss_phase`、`val/loss_class`
- `val/mean_phase_pearson`、`val/wrt_pearson`、`val/macro_f1`

### 4.4 checkpoint

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
- `phase_pearson_ES/MS/LS`、`phase_spearman_ES/MS/LS`
- `mean_phase_pearson`、`wrt_pearson`、`wrt_spearman`
- `macro_f1`、`balanced_accuracy`、`ev_l_auroc`

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

## 6. In silico mutagenesis

### 6.1 Saturation mutagenesis

对序列每个位置逐一突变为所有 3 种替代碱基，计算 delta WRT：

```bash
python scripts/virtual_knockout.py saturation \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --fasta /path/to/genome.fa \
    --chrom chr1 --start 1000000 --end 1032768 \
    --output outputs/interpretation/saturation.tsv
```

### 6.2 Motif scramble

将指定 motif 区域随机打乱，测试 motif 对 RT 的贡献：

```bash
python scripts/virtual_knockout.py motif \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --fasta /path/to/genome.fa \
    --chrom chr1 --start 1000000 --end 1032768 \
    --motif_start 4000 --motif_end 4020 \
    --mut_mode scramble \
    --output outputs/interpretation/cpp_motif.tsv
```

---

## 7. 基线模型

```python
import sys; sys.path.insert(0, '.')
from src.baselines import run_gc_baseline, run_kmer_baseline

gc_result = run_gc_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class)
kmer_result = run_kmer_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class, k=4)
```

---

## 8. 项目结构

```
src/
  data/
    data_utils.py       # GenomeSequence, one_hot_encode, TPM norm, WRT, RT class, labels
    dataset.py          # SpeciesConfig/Manifest, RepliSeqDataset
  models/
    model.py            # Basenji2Model, _Basenji2Trunk, _DilatedResidual, _SharedHead, MultiTaskLoss
  configs/
    transformer_wt.yaml # 训练超参数
  trainer.py            # DDP 训练循环 + TensorBoard
  eval.py               # 评估指标 + WT/cross-species 评估
  interpret.py          # saturation mutagenesis, motif mutation, IG
  baselines.py          # GC baseline, k-mer linear baseline

scripts/
  train.py              # 训练入口（支持 torchrun）
  evaluate.py           # 评估入口
  virtual_knockout.py   # in silico mutagenesis 入口

data/
  manifest.yaml         # 物种清单（路径 + 染色体划分）
  genomes/              # 参考基因组 FASTA（不入 git）
  labels/               # Repli-seq count TSV（不入 git）

outputs/                # 模型 checkpoint、评估结果（不入 git）
logs/                   # TensorBoard 日志（不入 git）
```

---

## 9. 模型框架

### 9.1 整体结构

```
DNA sequence [32 kb]
      │
 one_hot_encode
 [B, 4, L]  (A/C/G/T，N→全零)
      │
 Conv stem
 Conv1d(4→288, k=15, pool=2)          → L/2
      │
 Conv tower（6 层，每层 pool=2）
 Conv1d(288→339, k=5, pool=2)         → L/4
 Conv1d(339→399, k=5, pool=2)         → L/8
 Conv1d(399→470, k=5, pool=2)         → L/16
 Conv1d(470→554, k=5, pool=2)         → L/32
 Conv1d(554→652, k=5, pool=2)         → L/64
 Conv1d(652→768, k=5, pool=2)         → L/128 = 256
 每层: BN → GELU → Conv
      │
 Projection: Conv1d(768→384, k=1)
      │
 Dilated residual tower（11 层）
 dilation = round(1.5^i), i=0..10
 每层: [BN→GELU→Conv(k=3,dil)] + [BN→GELU→Conv(k=1)] + residual add
      │
 Cropping1D（两端各裁 16）
      │
 Bottleneck: BN→GELU→Conv1d(384→1536, k=1, dropout=0.05)
      │
 Global Average Pool → [B, 1536]
      │
 _SharedHead
 Linear(1536, 128) + GELU + Dropout
      ┌──┴──┐
 phase_out → Softplus → phase_pred [B, 3]  (ES/MS/LS log1p TPM)
 class_out → class_logits [B, 4]           (ES/MS/LS/NR)
```

### 9.2 Dilated Residual Block

每个 block 包含两个卷积，第一个带 dilation，第二个 k=1 做 projection，gamma 零初始化（InitZero trick）保证训练初期残差分支输出接近零：

```
x [B, C, L]
│
├─ BN → GELU → Conv(k=3, dilation=d)   ← 扩大感受野
│       │
│  BN → GELU → Conv(k=1, dropout=0.3)  ← gamma 零初始化
│
x = x + residual
```

11 层的 dilation 序列（rate_mult=1.5，取整）：1, 2, 2, 3, 5, 8, 11, 17, 25, 38, 57，最大感受野约 ±57 × 3 = ±171 个 128 bp bin，覆盖约 ±22 kb。

### 9.3 MultiTaskLoss

```
total = λ_phase × MSE(phase_pred, phase_labels)
      + λ_class × CrossEntropy(class_logits, rt_class)
```

默认 `λ_phase=1.0`，`λ_class=0.5`。`class_weights` 由训练集各类频率的平方根倒数计算，处理类别不平衡。早停以 `val_loss_total` 为准。

### 9.4 训练策略（对齐 Basenji2）

| 参数 | 值 |
|------|----|
| Optimizer | SGD + Nesterov momentum |
| lr | 0.15 |
| momentum | 0.99 |
| Scheduler | ReduceLROnPlateau（patience=8，factor=0.2） |
| Clip norm | 2.0 |
| Batch size | 8 |
| Early stopping patience | 16 |
| 数据增强 | reverse complement（p=0.5） |
