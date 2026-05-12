# Enformer-style Repli-seq RT Classifier

基于 DNA 序列的 Enformer 风格模型，对多物种基因组每个 128 bp bin 预测复制时序类别（ES / MS / LS / NR）。架构严格对齐 Enformer：one-hot 输入 → conv stem → conv tower（attention pooling）→ Transformer tower（TransformerXL 相对位置编码）→ 每物种独立分类头。

---

## 目录

1. [环境安装](#1-环境安装)
2. [数据准备](#2-数据准备)
3. [添加新物种](#3-添加新物种)
4. [配置文件](#4-配置文件)
5. [训练](#5-训练)
6. [评估](#6-评估)
7. [Checkpoint 迁移](#7-checkpoint-迁移)
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

### 2.2 Raw count 前处理

原始输入为 BED 格式的 bin × sample count 矩阵（128 bp 和 1024 bp 各一份），由 `scripts/preprocess_repliseq.py` 处理后生成 GFF3 标签文件。

#### 输入格式

无表头的 TSV，列顺序为 `chr  start  end  G1  ES_rep1  ES_rep2  ...`：

```
1   0       128     45   12   18   20   ...
1   128     256     312  89   95   101  ...
```

- 坐标为 0-based 半开区间（BED 标准）
- 列顺序需与脚本中 `SAMPLE_COLS` 定义一致
- 染色体名称需与 FASTA 中的 sequence name 一致

#### 标签分配流程

```
raw counts
    │
    ▼
CPM 归一化（全基因组 library size）
    │
    ▼
各 S 期窗口内 replicate 取平均
    │
    ▼
log2((S_avg + 1) / (G1 + 1))   ← 每个 bin 得到 ES / MS / LS 三个 score
    │
    ├── max(ES, MS, LS) ≤ 0  →  Non-replication (NR)
    └── max(ES, MS, LS) > 0  →  argmax → ES / MS / LS
```

**关键参数：**
- `NON_REP_THRESHOLD = 0.0`：log2 ratio ≤ 0 的 bin 判定为 NR（S 期无富集）
- CPM 归一化在全基因组范围内计算 library size，不按染色体分别计算

#### 运行前处理

编辑 `scripts/preprocess_repliseq.py` 中的路径和列定义，然后运行：

```bash
python scripts/preprocess_repliseq.py
```

输出到 `data/labels/`，每个物种 × 分辨率生成两个文件：
- `<species>_1024bp.gff3` — 用于训练（填入 manifest）
- `<species>_128bp.gff3` — 备用（更细粒度分析）
- `<species>_1024bp.tsv` / `<species>_128bp.tsv` — 含 label 和 score 列，供 QC 使用

#### 各物种列定义

**拟南芥（arabidopsis）**：G1 + ES/MS/LS 各 3 个 replicate，共 10 列

```python
ARA_SAMPLE_COLS = [
    "G1",
    "ES_rep1", "ES_rep2", "ES_rep3",
    "MS_rep1", "MS_rep2", "MS_rep3",
    "LS_rep1", "LS_rep2", "LS_rep3",
]
```

**玉米（zeamay/maize）**：G1 + ES/MS 各 3 个 replicate + LS 2 个 replicate，共 9 列

```python
ZEA_SAMPLE_COLS = [
    "G1",
    "ES_rep1", "ES_rep2", "ES_rep3",
    "MS_rep1", "MS_rep2", "MS_rep3",
    "LS_rep1", "LS_rep2",            # 无 LS_rep3
]
```

脚本自动处理缺失 replicate（取现有 replicate 的均值）。

#### 输出 GFF3 格式

```
##gff-version 3
1   assign_rt_class   peaks   1      1024   .   .   .   Name=ES;color=#2250F1;
1   assign_rt_class   peaks   1025   2048   .   .   .   Name=MS;color=#9B30FF;
1   assign_rt_class   peaks   2049   3072   .   .   .   Name=Non-replication;color=#B0B0B0;
```

- 坐标为 1-based 闭区间（GFF3 标准），由 BED 坐标自动转换
- `feature` 字段固定为 `peaks`，`Name=` 属性为 RT 类别

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
    gff3: /path/to/labels.gff3
    train_chroms: [chr01, chr02, chr03, chr05, chr06, chr07, chr09, chr10]
    val_chroms: [chr04, chr08]
    test_chroms: [chr11, chr12]
```

**重要：`species_id` 必须显式指定且全局唯一。** 不要依赖列表顺序——在中间插入新物种不会影响已有物种的 ID，checkpoint 路由不会错位。

当前已配置物种：

| 物种 | species_id | train | val | test |
|------|-----------|-------|-----|------|
| 水稻 (rice) | 0 | chr01–03, 05–07, 09–10 | chr04, chr08 | chr11, chr12 |
| 玉米 (maize) | 1 | chr1–4, 6–8 | chr5, chr9 | chr10 |
| 拟南芥 (arabidopsis) | 2 | chr1–3 | chr4 | chr5 |

---

## 4. 配置文件

### 4.1 manifest（物种清单）

`data/manifest.yaml` — 见第 3 节。

### 4.2 训练配置

`src/configs/transformer_wt.yaml` 包含所有超参数，关键字段：

```yaml
data:
  manifest: data/manifest.yaml
  input_window_length: 196608    # 196 kb 输入窗口，对齐 Enformer

augmentation:
  rc_prob: 0.5                   # reverse complement 概率
  shift_max: 1024                # 训练时随机平移最大 bp 数

model:
  bn_momentum: 0.1               # TF bn_momentum=0.9 → PyTorch 1-0.9=0.1

loss:
  class_weights: [1.43, 1.43, 1.43, 1.0]  # ES, MS, LS, NR

training:
  learning_rate: 0.0001
  warmup_steps: 2000
  batch_size: 4
  gradient_accumulation_steps: 8
  gradient_clip_norm: 0.2
  mixed_precision: true
  early_stopping_patience: 10
  max_epochs: 100
  seed: 42
```

---

## 5. 训练

### 5.1 多物种训练策略

多物种训练的核心思路：**trunk 参数在所有物种间共享，每个物种有独立的输出 head**。训练时所有物种的样本混合在同一个 batch 中，按 `species_id` 分组路由到各自的 head。

#### 采样均衡

单卡训练使用 `SpeciesBalancedSampler`，保证每个 batch 中各物种样本数相等：

```
batch_size = 4,  3 个物种
→ 每物种 1 个样本/batch（batch_size // num_species）
→ 样本少的物种自动循环过采样
→ 训练轮数由样本最多的物种决定
```

这样 trunk 在每个 step 都能同时接收来自所有物种的梯度，避免某一物种主导参数更新。

#### 损失计算

每个 batch 按物种分组，分别计算 cross-entropy，再按样本数加权平均：

```
batch B = {rice: n_r 个样本, maize: n_m 个样本, arabidopsis: n_a 个样本}

loss_rice  = CE(logits_rice,  labels_rice)   # 对 n_r 个样本的 896 个 bin 求均值
loss_maize = CE(logits_maize, labels_maize)
loss_ara   = CE(logits_ara,   labels_ara)

loss_total = (n_r × loss_rice + n_m × loss_maize + n_a × loss_ara) / B
```

加权平均而非简单平均，确保 loss 的量纲与 batch size 无关，梯度累积时数值稳定。

CE 使用类别权重 `[1.43, 1.43, 1.43, 1.0]`（ES / MS / LS / NR），补偿 NR 类别在基因组中占比偏高的问题。

#### 梯度流向

```
loss_total.backward()
    │
    ├── ∂loss/∂head_rice    → 只更新 rice head
    ├── ∂loss/∂head_maize   → 只更新 maize head
    ├── ∂loss/∂head_ara     → 只更新 arabidopsis head
    └── ∂loss/∂trunk        → 所有物种的梯度叠加，共同更新 trunk
```

trunk 接收来自所有物种的梯度信号，学习跨物种通用的序列特征表示。

### 5.2 单卡训练

```bash
python scripts/train.py --config src/configs/transformer_wt.yaml
```

### 5.3 多卡 DDP 训练

```bash
torchrun --nproc_per_node=4 scripts/train.py --config src/configs/transformer_wt.yaml
```

DDP 模式使用 `DistributedSampler`（`SpeciesBalancedSampler` 暂不支持多卡，各 rank 的 batch 物种分布可能不均衡）。

### 5.4 断点续训

```bash
python scripts/train.py \
    --config src/configs/transformer_wt.yaml \
    --resume outputs/checkpoints/best_model.pt
```

### 5.5 训练过程监控

```bash
tensorboard --logdir logs/
```

TensorBoard 记录：
- `train/loss_total`、`train/loss_rt`（epoch 平均）
- `train/lr`
- `val/loss_total`、`val/loss_rt`
- `val/macro_f1`、`val/overall_acc`
- `val/acc_ES`、`val/acc_MS`、`val/acc_LS`、`val/acc_NR`

### 5.6 Checkpoint

最优模型保存在 `outputs/checkpoints/best_model.pt`，由 `val_loss_total`（越低越好）决定。Checkpoint 包含所有物种的 head 权重，可直接用于新物种的 fine-tuning（新物种 head 随机初始化，trunk 从 checkpoint 加载）。

---

## 6. 评估

### 6.1 Per-species 测试集评估

对每个物种使用其自己的 head 在测试集上评估：

```bash
python scripts/evaluate.py \
    --mode wt \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/wt_test_metrics.tsv`，每行一个物种，包含：
- `overall_acc`、`macro_f1`
- `acc_ES`、`acc_MS`、`acc_LS`、`acc_NR`
- `cm_row_ES`、`cm_row_MS`、`cm_row_LS`、`cm_row_NR`（混淆矩阵各行）

### 6.2 N×N 跨物种评估

用所有物种的 head 分别预测所有物种的测试集，输出 N×N 矩阵：

```bash
python scripts/evaluate.py \
    --mode cross_species \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/cross_species_metrics.tsv`，包含列 `src_head`（head 来源物种）和 `tgt_species`（测试数据物种）。对角线为 per-species 结果，非对角线为跨物种迁移能力。

---

## 7. Checkpoint 迁移

如果你有旧版单 head 的 checkpoint（`head.weight` / `head.bias`），需要迁移到新的 ModuleDict 格式：

```bash
python scripts/migrate_checkpoint.py \
    old_checkpoint.pt \
    new_checkpoint.pt \
    --species rice
```

迁移后的 checkpoint 可直接用于多物种训练（rice head 已初始化，其他物种 head 随机初始化）。

---

## 8. 项目结构

```
src/
  data/
    data_utils.py         # GenomeSequence, one_hot_encode, load_labels, load_labels_indexed
    dataset.py            # SpeciesConfig, load_manifest, RepliSeqDataset, SpeciesBalancedSampler
  models/
    model.py              # Basenji2Model（ModuleDict heads）, RTClassLoss
  configs/
    transformer_wt.yaml   # 训练超参数
  trainer.py              # DDP 训练循环，_forward_multi_species，TensorBoard
  eval.py                 # evaluate_predictions, eval_wt, eval_cross_species

scripts/
  train.py                # 训练入口（支持 torchrun）
  evaluate.py             # 评估入口（--mode wt | cross_species）
  migrate_checkpoint.py   # 单 head → ModuleDict 迁移工具
  preprocess_repliseq.py  # Repli-seq 数据预处理

data/
  manifest.yaml           # 物种清单（路径 + 染色体划分 + species_id）
  genomes/                # 参考基因组 FASTA（不入 git）
  labels/                 # GFF3 标签文件（不入 git）

outputs/                  # 模型 checkpoint、评估结果（不入 git）
logs/                     # TensorBoard 日志（不入 git）
tests/                    # 单元测试
docs/
  superpowers/
    specs/                # 设计文档
    plans/                # 实现计划
```

---

## 9. 模型框架

### 9.1 整体结构

```
DNA sequence [196 kb]
[B, 4, 196608]
      │
 Conv stem: Conv1d(4→768, k=15) + Residual(ConvBlock) + AttnPool/2
[B, 768, 98304]
      │
 Conv tower（6 层，每层 AttnPool/2）
 filters: 768→768→896→1024→1152→1280→1536
[B, 1536, 1536]  @ 128 bp/token
      │
 Transformer tower（11 层，d_model=1536, n_heads=8, dim_key=64）
 TransformerXL 相对位置编码（exponential + central_mask + gamma basis）
[B, 1536, 1536]
      │
 Crop 320 each side → [B, 896, 1536]
      │
 Final pointwise: BN → GELU → Conv1d(1536→3072) → Dropout(0.05) → GELU
[B, 896, 3072]
      │
 Per-species head: Linear(3072→4)   ← 每物种独立，共享 trunk
[B, 896, 4]  → ES / MS / LS / NR per 128 bp bin
```

### 9.2 多物种设计

模型使用 `nn.ModuleDict` 管理每个物种的独立输出头：

```python
self._heads = nn.ModuleDict({
    "rice":        nn.Linear(3072, 4),
    "maize":       nn.Linear(3072, 4),
    "arabidopsis": nn.Linear(3072, 4),
})
```

- **Trunk 完全共享**：所有物种共用同一套 conv + transformer 参数
- **Head 完全独立**：每个物种有自己的 `Linear(3072→4)`，可学习物种特异性偏差
- **Forward 路由**：`model(one_hot, head="rice")` 指定使用哪个 head
- **新增物种**：只需在 `manifest.yaml` 加条目，模型初始化时自动创建对应 head

### 9.3 训练策略

混合 batch 训练：每个 batch 包含所有物种的样本（`SpeciesBalancedSampler` 保证均衡），按 `species_id` 分组路由到对应 head，加权平均 loss：

```
loss = Σ (n_sp / B) × CE(logits_sp, labels_sp)
```

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

### 9.4 采样策略

滑动窗口，stride = 输出区域大小（零重叠，对齐 Enformer）：

```
stride = 896 × 128 = 114688 bp
每个 bin 只出现在一个训练样本中
```

标签坐标对齐：第一个输出 bin 对应 `win_start + 320 × 128 = win_start + 40960 bp`。
