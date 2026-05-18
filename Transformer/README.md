# Repliformer — Repli-seq RT 信号回归模型

基于 DNA 序列的 Enformer 风格模型，对多物种基因组每个 128 bp bin 回归预测复制时序（RT）连续信号（G1 / ES / MS / LS 四通道）。架构严格对齐 Enformer：one-hot 输入 → conv stem → conv tower（attention pooling）→ Transformer tower（TransformerXL 相对位置编码）→ 每物种独立回归头。所有架构参数通过 `RepliformerConfig` 集中管理，可在配置文件中直接修改。

---

## 目录

1. [环境安装](#1-环境安装)
2. [数据准备](#2-数据准备)
3. [添加新物种](#3-添加新物种)
4. [配置文件](#4-配置文件)
5. [训练](#5-训练)
6. [评估](#6-评估)
7. [项目结构](#7-项目结构)
8. [模型框架](#8-模型框架)

---

## 1. 环境安装

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install pyfaidx pandas numpy scipy scikit-learn pyyaml tqdm tensorboard transformers
```

Python 3.10+，PyTorch 2.x。

---

## 2. 数据准备

### 2.1 参考基因组

将各物种的参考基因组 FASTA 文件放到对应目录，并建立 pyfaidx 索引：

```bash
python -c "from pyfaidx import Fasta; Fasta('/path/to/genome.fa')"
```

### 2.2 信号前处理

原始输入为 BED 格式的 bin × sample count 矩阵（128 bp），由 `data/preprocess_signals.py` 处理后生成模型所需的 TSV 信号文件。

#### 输入格式

无表头的 TSV，列顺序为 `chr  start  end  sample1  sample2  ...`（各物种列定义见脚本中 `SPECIES` 字典）。

#### 处理流程

```
raw counts
    │
    ▼
CPM 归一化（counts per million，消除 library size 差异）
    │
    ▼
每个 replicate 单独 soft clip（x > 384 → 383 + sqrt(x - 383)）
    │
    ▼
各期 replicate 取平均 → G1, ES, MS, LS
    │
    ▼
log1p 变换（压缩动态范围，保持非负）
    │
    ▼
输出 TSV: chrom  start  end  G1  ES  MS  LS
```

各物种 head 独立学习各自的信号尺度，不需要跨物种归一化。

#### 运行前处理

```bash
# 处理所有物种
python data/preprocess_signals.py

# 处理单个物种
python data/preprocess_signals.py rice
python data/preprocess_signals.py maize arabidopsis
```

输出到 `data/labels/`：
- `rice_128bp_rt_signals.tsv`
- `zeamay_128bp_rt_signals.tsv`
- `arabidopsis_128bp_rt_signals.tsv`

#### 各物种列定义

**水稻（rice）**：G1 × 2 replicate，ES × 2，MS × 1，LS × 1

**玉米（maize）**：G1 × 1，ES × 3，MS × 3，LS × 2

**拟南芥（arabidopsis）**：G1 × 1，ES × 3，MS × 3，LS × 3

#### 输出 TSV 格式

```
chrom   start   end     G1        ES        MS        LS
chr01   0       128     0.000000  0.000000  0.000000  0.000000
chr01   896     1024    0.279134  0.312456  0.198234  0.267891
```

- 坐标为 0-based 半开区间（BED 标准）
- 四列信号为 `log1p(soft_clip(CPM))`，非负 float32

### 2.3 数据划分说明

划分以**完整染色体**为单位，禁止随机划分（相邻 bin 序列高度相关，随机划分会造成 genomic leakage）。

---

## 3. 添加新物种

**只需编辑 `data/manifest.yaml`，无需改动任何代码。**

每个物种条目包含以下字段：

```yaml
species:
  - name: rice           # 物种名（唯一标识符，用作 head key）
    species_id: 0        # 整数 ID，必须唯一，不依赖列表顺序
    fasta: /path/to/genome.fasta
    tsv: /path/to/rt_signals.tsv
    train_chroms: [chr01, chr02, chr03, chr05, chr06, chr07, chr09, chr10]
    val_chroms: [chr04, chr08]
    test_chroms: [chr11, chr12]
```

**重要：`species_id` 必须显式指定且全局唯一。** 不要依赖列表顺序——在中间插入新物种不会影响已有物种的 ID，checkpoint 路由不会错位。

当前已配置物种：

| 物种 | species_id | train | val | test |
|------|-----------|-------|-----|------|
| 水稻 (rice) | 0 | chr01, 03–09 | chr02, chr10 | chr11, chr12 |
| 玉米 (maize) | 1 | 1, 4–9 | 2, 3 | 10 |
| 拟南芥 (arabidopsis) | 2 | 1, 3, 4 | 2 | 5 |

---

## 4. 配置文件

### 4.1 manifest（物种清单）

`data/manifest.yaml` — 见第 3 节。

### 4.2 训练配置

`src/configs/transformer_wt.yaml` 包含所有超参数。`model:` 节的每个字段都直接对应 `RepliformerConfig` 的参数，修改后重启训练即生效，无需改动代码。

```yaml
data:
  manifest: data/manifest.yaml
  input_window_length: 196608    # 196 kb 输入窗口，对齐 Enformer

augmentation:
  rc_prob: 0.5                   # reverse complement 概率
  shift_max: 1024                # 训练时随机平移最大 bp 数

model:
  # Conv trunk
  pool_type: attn                # "attn" | "avg_late"
  stem_channels: 256
  stem_kernel_size: 15
  tower_filter_list: [256, 256, 320, 384, 384, 512, 512]
  bn_momentum: 0.1               # TF bn_momentum=0.9 → PyTorch 1-0.9=0.1
  # Transformer tower
  d_model: 384
  n_heads: 8
  dim_key: 64
  n_layers: 2
  dropout: 0.1
  attn_dropout: 0.05
  pos_dropout: 0.01
  # Final pointwise
  pointwise_channels: 1024
  final_dropout: 0.1
  # Sequence geometry
  crop_size: 320                 # 每侧裁剪的 bin 数（输出 896 bins）
  n_output_bins: 4               # 信号通道数：G1, ES, MS, LS

loss:
  type: mse

training:
  learning_rate: 0.0005
  weight_decay: 0.01
  warmup_steps: 1000
  batch_size: 4
  gradient_accumulation_steps: 4
  gradient_clip_norm: 1.0
  mixed_precision: true
  early_stopping_patience: 32
  max_epochs: 200
  seed: 42
```

---

## 5. 训练

### 5.1 多物种训练策略

**trunk 参数在所有物种间共享，每个物种有独立的输出 head。** 训练时各物种分别跑完整的 epoch，按物种轮流迭代，每个 step 只处理一个物种的 batch。

#### 损失函数

MSE loss，对每个 bin 的四通道信号回归：

```
loss = MSE(pred_signals, true_signals)   # NaN bin 自动 mask
```

模型输出经 `softplus` 激活保证非负，与 `log1p(CPM)` target 空间匹配。

#### WRT

WRT（Weighted Replication Timing）由预测的 ES/MS/LS 信号推导：

```
WRT = (0.5 × MS + LS) / (ES + MS + LS + ε)
```

在验证时自动计算并报告 WRT Pearson，无需额外预测头。

#### 梯度流向

```
loss.backward()
    ├── ∂loss/∂head_{species}  → 只更新当前物种的 head
    └── ∂loss/∂trunk           → 更新共享 trunk
```

### 5.2 单卡训练

```bash
python scripts/train.py --config src/configs/transformer_wt.yaml
```

### 5.3 多卡 DDP 训练

```bash
torchrun --nproc_per_node=4 scripts/train.py --config src/configs/transformer_wt.yaml
```

### 5.4 断点续训

```bash
python scripts/train.py \
    --config src/configs/transformer_wt.yaml \
    --resume outputs/basenji2_wt/checkpoints/best_model.pt
```

### 5.5 训练过程监控

```bash
tensorboard --logdir logs/
```

TensorBoard 记录：
- `train/loss_total`、`train/loss_rt`（epoch 平均）
- `train/lr`
- `val/loss_total`、`val/loss_rt`
- `val/pearson_mean`（所有物种平均 Pearson）
- `val/{species}/pearson_{G1,ES,MS,LS,WRT}`（每物种每通道 + WRT Pearson）

日志示例：

```
epoch=76 time=817s train_loss=0.0016 val_loss=0.0032 val_pearson=0.6210
  per-species Pearson: rice=0.618 maize=0.624 arabidopsis=0.623
    rice: ES=0.631 MS=0.625 LS=0.612 G1=0.604 WRT=0.651
    maize: ES=0.638 MS=0.621 LS=0.618 G1=0.619 WRT=0.643
    arabidopsis: ES=0.629 MS=0.618 LS=0.615 G1=0.630 WRT=0.638
```

### 5.6 Checkpoint

最优模型保存在 `outputs/basenji2_wt/checkpoints/best_model.pt`，由 `val_loss_total`（越低越好）决定。Checkpoint 包含所有物种的 head 权重，可直接用于新物种的 fine-tuning（新物种 head 随机初始化，trunk 从 checkpoint 加载）。

---

## 6. 评估

### 6.1 Per-species 测试集评估

对每个物种使用其自己的 head 在测试集上评估：

```bash
python scripts/evaluate.py \
    --mode wt \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/wt_test_metrics.tsv`，每行一个物种，包含：
- `mean_pearson`
- `pearson_G1`、`pearson_ES`、`pearson_MS`、`pearson_LS`、`pearson_WRT`
- `mse`、`mae`

### 6.2 N×N 跨物种评估

用所有物种的 head 分别预测所有物种的测试集，输出 N×N 矩阵：

```bash
python scripts/evaluate.py \
    --mode cross_species \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/cross_species_metrics.tsv`，包含列 `src_head`（head 来源物种）和 `tgt_species`（测试数据物种）。对角线为 per-species 结果，非对角线为跨物种迁移能力。

---

## 7. 项目结构

```
src/
  data/
    data_utils.py         # GenomeSequence, one_hot_encode, load_signals, load_signals_indexed
    dataset.py            # SpeciesConfig, load_manifest, RepliSeqDataset
  models/
    model.py              # RepliformerModel（ModuleDict heads）, RTSignalLoss
    config_model.py       # RepliformerConfig（所有架构参数）
  configs/
    transformer_wt.yaml   # 训练超参数
  trainer.py              # 训练循环，_validate，TensorBoard
  eval.py                 # evaluate_predictions（含 WRT）, eval_wt, eval_cross_species

scripts/
  train.py                # 训练入口（支持 torchrun）
  evaluate.py             # 评估入口（--mode wt | cross_species）
  migrate_checkpoint.py   # 单 head → ModuleDict 迁移工具

data/
  manifest.yaml           # 物种清单（路径 + 染色体划分 + species_id）
  preprocess_signals.py   # Repli-seq 信号前处理（raw counts → log1p CPM TSV）
  genomes/                # 参考基因组 FASTA（不入 git）
  labels/                 # TSV 信号文件（不入 git）

outputs/                  # 模型 checkpoint、评估结果（不入 git）
logs/                     # TensorBoard 日志（不入 git）
tests/                    # 单元测试
```

---

## 8. 模型框架

### 8.1 整体结构

```
DNA sequence [196 kb]
[B, 4, 196608]
      │
 Conv stem: Conv1d(4→256, k=15) + Residual(ConvBlock) + AttnPool/2
[B, 256, 98304]
      │
 Conv tower（6 层，每层 AttnPool/2）
 filters: 256→256→320→384→384→512→512
[B, 512, 1536]  @ 128 bp/token
      │
 trunk_proj: Linear(512→384)
      │
 Transformer tower（2 层，d_model=384, n_heads=8, dim_key=64）
 TransformerXL 相对位置编码（exponential + central_mask + gamma basis）
[B, 1536, 384]
      │
 Crop 320 each side → [B, 896, 384]
      │
 Final pointwise: BN → GELU → Conv1d(384→1024) → Dropout(0.1) → GELU
[B, 896, 1024]
      │
 Per-species head: Linear(1024→4) → softplus
[B, 896, 4]  → G1 / ES / MS / LS，每 128 bp bin 一组
      │
 WRT（推导，不预测）= (0.5×MS + LS) / (ES + MS + LS + ε)
```

### 8.2 RepliformerConfig

所有架构超参数集中在 `src/models/config_model.py` 的 `RepliformerConfig` 中，继承自 `transformers.PretrainedConfig`。训练/评估时从 yaml 的 `model:` 节自动构建：

```python
model_cfg = RepliformerConfig(**cfg.get("model", {}))
model = RepliformerModel(species_configs, model_cfg=model_cfg)
```

未在 yaml 中指定的字段使用 `RepliformerConfig` 的默认值。

### 8.3 多物种设计

模型使用 `nn.ModuleDict` 管理每个物种的独立输出头：

```python
self._heads = nn.ModuleDict({
    "rice":        nn.Linear(cfg.pointwise_channels, cfg.n_output_bins),
    "maize":       nn.Linear(cfg.pointwise_channels, cfg.n_output_bins),
    "arabidopsis": nn.Linear(cfg.pointwise_channels, cfg.n_output_bins),
})
```

- **Trunk 完全共享**：所有物种共用同一套 conv + transformer 参数
- **Head 完全独立**：每个物种有自己的 `Linear(1024→4)`，学习物种特异性尺度和偏差
- **Forward 路由**：`model(one_hot, head="rice")` 指定使用哪个 head
- **新增物种**：只需在 `manifest.yaml` 加条目，模型初始化时自动创建对应 head

### 8.4 训练超参数

| 参数 | 值 |
|------|----|
| Optimizer | AdamW |
| lr | 5e-4 |
| weight_decay | 0.01 |
| Scheduler | warmup (1000 steps) + cosine decay |
| Clip norm | 1.0 |
| 有效 batch size | 16（batch 4 × 梯度累积 4） |
| 数据增强 | reverse complement (p=0.5) + random shift (±1024 bp) |
| Loss | MSE |
| 输出激活 | softplus（非负，匹配 log1p CPM target） |
| Target 变换 | log1p(soft_clip(CPM))，预处理时完成 |

### 8.5 采样策略

滑动窗口，stride = 输出区域大小（零重叠）：

```
stride = (196608 - 2 × 320 × 128) = 114688 bp
每个 bin 只出现在一个训练样本中
```

标签坐标对齐：第一个输出 bin 对应 `win_start + 320 × 128 = win_start + 40960 bp`。
