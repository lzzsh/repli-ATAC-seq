# Enformer-style Repli-seq RT Classifier

基于 DNA 序列的 Enformer 风格模型，对水稻基因组每个 128 bp bin 预测复制时序类别（ES / MS / LS / NR）。架构对齐 Enformer：one-hot 输入 → conv stem → conv tower → bottleneck → Transformer tower → 逐 bin 分类头。

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
  input_window_length: 196608  # 196 kb 输入窗口，对齐 Enformer

augmentation:
  rc_prob: 0.5                 # reverse complement 概率
  shift_max: 1024              # 训练时随机平移最大 bp 数

model:
  bn_momentum: 0.1             # TF bn_momentum=0.9 → PyTorch 1-0.9=0.1

loss:
  class_weights: [1.43, 1.43, 1.43, 1.0]  # ES, MS, LS, NR

training:
  learning_rate: 0.0001        # Adam lr
  warmup_steps: 2000           # 线性 warmup 步数
  batch_size: 4                # 等效 batch = 4 × 8 = 32（梯度累积）
  gradient_clip_norm: 0.2
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
DNA sequence [196 kb]
[B, 4, 196608]
      │
 Conv stem: Conv1d(4→288, k=15) + AttnPool/2
[B, 288, 98304]
      │
 Conv tower（6 层，每层 AttnPool/2）
 filters: 288→339→399→470→554→652→768  (×1.1776 递增, k=5)
[B, 768, 1536]  @ 128 bp/token
      │
 Crop 320 each side
[B, 768, 896]
      │
 Bottleneck: Conv1d(768→1536, k=1, dropout=0.05)
[B, 1536, 896]
      │
 Transformer tower（11 层，d_model=1536, n_heads=8, ffn=3072, dropout=0.4）
 + 相对位置偏置（Enformer-style learnable bias）
[B, 896, 1536]
      │
 LayerNorm + Linear(1536→4)
[B, 896, 4]  → ES / MS / LS / NR per 128 bp bin
```

### 7.2 采样策略

滑动窗口，stride = 输出区域大小（零重叠，对齐 Enformer）：

```
stride = 896 × 128 = 114688 bp
每个 bin 只出现在一个训练样本中
```

标签坐标对齐：第一个输出 bin 对应 `win_start + 320 × 128 = win_start + 40960 bp`。

### 7.3 训练策略（对齐 Enformer）

| 参数 | 值 |
|------|----|
| Optimizer | Adam |
| lr | 1e-4 |
| betas | (0.9, 0.999) |
| Scheduler | warmup (2000 steps) + cosine decay |
| Clip norm | 0.2 |
| 有效 batch size | 32（4 × 梯度累积 8） |
| 数据增强 | reverse complement (p=0.5) + random shift (±1024 bp) |
| Loss | cross-entropy + class weights [1.43, 1.43, 1.43, 1.0] |
