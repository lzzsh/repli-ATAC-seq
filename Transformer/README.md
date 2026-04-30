# Transformer Repli-seq

基于 DNA 序列的多物种 Transformer 模型，预测 Repli-seq ES/MS/LS 信号和 RT class，支持 in silico mutagenesis 解释 CPP motif 对复制时序的贡献。架构对齐 Enformer/Basenji2：one-hot 输入 → 6 层 conv tower（64x 下采样）→ Transformer encoder。

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
mkdir -p data/genomes/rice data/genomes/maize data/genomes/arabidopsis

# 放入 FASTA 文件后建索引（pyfaidx 会自动建，也可手动）
python -c "from pyfaidx import Fasta; Fasta('data/genomes/rice/IRGSP-1.0.fa')"
python -c "from pyfaidx import Fasta; Fasta('data/genomes/maize/Zm-B73-REFERENCE-NAM-5.0.fa')"
python -c "from pyfaidx import Fasta; Fasta('data/genomes/arabidopsis/TAIR10.fa')"
```

### 2.2 Repli-seq count 文件

每个物种需要一个 TSV 文件，每行对应一个 1 kb bin，格式如下：

```
chrom   start   end     ES_count    MS_count    LS_count
chr1    0       1000    142         87          23
chr1    1000    2000    98          201         56
...
```

- `chrom`：染色体名称，需与 FASTA 中的 sequence name 一致
- `start` / `end`：bin 的基因组坐标（0-based，半开区间）
- `ES_count` / `MS_count` / `LS_count`：各 phase 的 raw read count

将文件放到：

```
data/labels/rice_wt_counts.tsv
data/labels/maize_wt_counts.tsv
data/labels/arabidopsis_wt_counts.tsv
```

> **归一化由代码自动完成**：TPM 归一化 → log1p 变换 → WRT 计算 → RT class 分配（E/M/L）。EM/ML 边界 bins 会被自动过滤。

### 2.3 验证数据格式

```bash
python -c "
import sys; sys.path.insert(0, '.')
from src.data.data_utils import load_labels
df = load_labels('data/labels/rice_wt_counts.tsv', 'rice')
print(df.shape)
print(df[['chrom','start','end','WRT','RT_class']].head())
"
```

### 2.4 数据划分说明

划分以**完整染色体**为单位，禁止随机划分（相邻 bins 序列高度相关，随机划分会造成 genomic leakage）。

默认划分（见 `data/manifest.yaml`）：

| 物种 | train | val | test |
|------|-------|-----|------|
| 水稻 | chr1–chr8 | chr9–chr10 | chr11–chr12 |
| 玉米 | chr1–chr7 | chr8 | chr9–chr10 |
| 拟南芥 | chr1–chr3 | chr4 | chr5 |

---

## 3. 配置文件

### 3.1 manifest（物种清单）

`data/manifest.yaml` 定义每个物种的数据路径和染色体划分：

```yaml
species:
  - name: rice
    fasta: data/genomes/rice/IRGSP-1.0.fa
    count_tsv: data/labels/rice_wt_counts.tsv
    train_chroms: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8]
    val_chroms: [chr9, chr10]
    test_chroms: [chr11, chr12]
  - name: maize
    ...
```

### 3.2 训练配置

`src/configs/transformer_wt.yaml` 包含所有超参数，关键字段：

```yaml
data:
  input_window_length: 131072   # 131 kb 输入窗口，对齐 Basenji2

model:
  d_model: 256        # 模型维度
  n_layers: 6         # Transformer 层数
  n_heads: 8          # 注意力头数
  dim_feedforward: 1024

training:
  learning_rate: 0.0001
  batch_size: 16
  max_epochs: 100
  early_stopping_patience: 10   # val_loss 连续 10 个 epoch 不下降则停止
  mixed_precision: true         # 开启 AMP，节省显存
```

---

## 4. 训练

### 4.1 单卡训练

```bash
python scripts/train.py --config src/configs/transformer_wt.yaml
```

### 4.2 多卡 DDP 训练（推荐）

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
- `val/mean_phase_pearson`、`val/wrt_spearman`、`val/macro_f1`

### 4.4 checkpoint

最优模型保存在 `outputs/transformer_wt/checkpoints/best_model.pt`，由 `val_loss_total`（越低越好）决定。

---

## 5. 评估

### 5.1 WT 测试集评估（各物种分别报告）

```bash
python scripts/evaluate.py \
    --mode wt \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --output outputs/metrics
```

输出 `outputs/metrics/wt_test_metrics.tsv`，包含每个物种的：
- `phase_pearson_ES/MS/LS`、`phase_spearman_ES/MS/LS`
- `mean_phase_pearson`、`wrt_spearman`
- `macro_f1`、`balanced_accuracy`、`ev_l_auroc`

### 5.2 Leave-one-species-out 评估

```bash
python scripts/evaluate.py \
    --mode cross_species \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --held_out arabidopsis \
    --output outputs/metrics
```

输出 `outputs/metrics/cross_species_metrics.tsv`。

---

## 6. In silico mutagenesis

需要提供目标区域的基因组坐标和对应物种的 FASTA 文件。

### 6.1 Saturation mutagenesis

对序列每个位置逐一突变为所有 3 种替代碱基，计算 delta WRT：

```bash
python scripts/virtual_knockout.py saturation \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --fasta data/genomes/rice/IRGSP-1.0.fa \
    --species_id 0 --n_species 3 \
    --chrom chr1 --start 1000000 --end 1008192 \
    --output outputs/interpretation/saturation.tsv
```

输出：`saturation.tsv`（每个突变的 delta_ES/MS/LS/WRT）和 `saturation_importance.npy`（每个位置的最大 |delta_WRT|）。

### 6.2 Motif scramble

将指定 motif 区域随机打乱（保留碱基组成），测试 motif 对 RT 的贡献：

```bash
python scripts/virtual_knockout.py motif \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --fasta data/genomes/rice/IRGSP-1.0.fa \
    --species_id 0 --n_species 3 \
    --chrom chr1 --start 1000000 --end 1008192 \
    --motif_start 4000 --motif_end 4020 \
    --mut_mode scramble \
    --output outputs/interpretation/cpp_motif.tsv
```

`--mut_mode` 可选 `scramble`（随机打乱）或 `gc_matched_random`（GC 含量匹配的随机序列）。

### 6.3 Flank control

同时测试 motif 本身、flanking 区域、以及 motif+flank 联合突变的效果，用于区分 motif 贡献和局部序列背景：

```bash
python scripts/virtual_knockout.py flank \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --fasta data/genomes/rice/IRGSP-1.0.fa \
    --species_id 0 --n_species 3 \
    --chrom chr1 --start 1000000 --end 1008192 \
    --motif_start 4000 --motif_end 4020 \
    --output outputs/interpretation/flank_control.tsv
```

### 6.4 Integrated Gradients

计算每个碱基位置对 WRT 预测的归因分数：

```bash
python scripts/virtual_knockout.py ig \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/transformer_wt/checkpoints/best_model.pt \
    --fasta data/genomes/rice/IRGSP-1.0.fa \
    --species_id 0 --n_species 3 \
    --chrom chr1 --start 1000000 --end 1008192 \
    --output outputs/interpretation/ig.npy
```

---

## 7. 基线模型

用于对比，不需要 GPU：

```python
import sys; sys.path.insert(0, '.')
from src.baselines import run_gc_baseline, run_kmer_baseline

# train_seqs / test_seqs: list of DNA strings
# train_wrt / test_wrt: np.ndarray of WRT values
# train_class / test_class: np.ndarray of int (0=E, 1=M, 2=L)

gc_result = run_gc_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class)
kmer_result = run_kmer_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class, k=4)
print(gc_result)
print(kmer_result)
```

---

## 8. 项目结构

```
src/
  data/
    data_utils.py       # GenomeSequence, one_hot_encode, TPM norm, WRT, RT class, labels
    dataset.py          # SpeciesConfig/Manifest, RepliSeqDataset, SpeciesBalancedSampler
  models/
    model.py            # _ConvStem, RoPE encoder, FiLM, AttentionPooling, _SharedHead, DNATransformer, MultiTaskLoss
  configs/
    transformer_wt.yaml # 训练超参数
  trainer.py            # DDP 训练循环 + TensorBoard
  eval.py               # 评估指标 + WT/cross-species 评估
  interpret.py          # saturation mutagenesis, motif mutation, flank control, IG
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

## species_id 对照表

`species_id` 由 `data/manifest.yaml` 中的顺序决定（0-indexed）：

| species_id | 物种 |
|-----------|------|
| 0 | rice |
| 1 | maize |
| 2 | arabidopsis |

新增物种时，在 `manifest.yaml` 末尾追加，`n_species` 参数相应增加，**已有 species_id 不变**。

---

## 9. 模型框架

### 9.1 整体结构

```
DNA sequence [131 kb]
      │
 one_hot_encode
 [B, 4, L]  (A/C/G/T，N→全零)
      │
 _ConvStem（6层，对齐 Enformer/Basenji2）
 Conv1d(k=15, stride=2, 64 filters)   → L/2
 Conv1d(k=5,  stride=2, 96 filters)   → L/4
 Conv1d(k=5,  stride=2, 128 filters)  → L/8
 Conv1d(k=5,  stride=2, 160 filters)  → L/16
 Conv1d(k=5,  stride=2, 192 filters)  → L/32
 Conv1d(k=5,  stride=2, 256 filters)  → L/64 ≈ 1024
 每层: Conv1d → LayerNorm → GELU → Dropout
      │                        │
      │                   species_emb [B, species_emb_dim]
      │                        │
      └──────────┬─────────────┘
                 │         (species_emb 传入每一层)
       ┌─────────┴─────────┐
       │   × n_layers       │
       │  _EncoderLayer     │
       │  pre-norm + RoPE   │
       │  + FiLM            │
       └─────────┬─────────┘
                 │
            LayerNorm
                 │
         AttentionPooling
            [B, d_model]
                 │
         _SharedHead
         Linear(d_model, hidden) + GELU + Dropout
              ┌──┴──┐
        phase_out   → Softplus → phase_pred [B, 3]  (ES/MS/LS)
                              ↓
                    WRT = (0.5×MS + LS) / (ES+MS+LS)
                              ↓
                    class_logits [B, 3]  (E/M/L，由 WRT 阈值推算)
```

输入为 one-hot 编码（对齐 Enformer/Basenji2），conv tower 6 层 64x 下采样后序列长度约 1024，进入 Transformer encoder。物种信息通过 FiLM 在每一层主动调制序列特征。

### 9.2 _EncoderLayer（pre-norm + RoPE + FiLM）

```
x [B, L, d_model]
│
├─ LayerNorm → FiLM(species_emb)     ← 物种调制 pre-attention 特征
│       │
│   Q/K/V proj（各 [B, L, d_model]）
│       │
│   reshape → [B, L, n_heads, head_dim]
│       │
│   _apply_rope(Q), _apply_rope(K)   ← RoPE 旋转位置编码
│       │
│   transpose → [B, n_heads, L, head_dim]
│       │
│   scores = Q @ Kᵀ * scale          ← scaled dot-product
│   masked_fill(key_padding_mask)    ← 可选 padding mask
│   softmax → attn [B, n_heads, L, L]
│       │
│   out = attn @ V → reshape → out_proj
│
x = x + out                          ← 残差
│
├─ LayerNorm → FiLM(species_emb)     ← 物种调制 pre-FFN 特征
│       │
│   FFN: Linear → GELU → Dropout → Linear → Dropout
│
x = x + FFN(x)                       ← 残差
```

### 9.3 RoPE 实现

```
_precompute_freqs_cis(head_dim, max_seq_len)
  → freqs_cis [max_seq_len, head_dim//2]  (complex64)
  → 注册为 non-persistent buffer，随模型自动迁移设备

forward 时：freqs_cis = self.freqs_cis[:L]  ← 按实际序列长度切片

_apply_rope(x, freqs_cis):
  x [B, L, n_heads, head_dim]
  → view_as_complex → [B, L, n_heads, head_dim//2]  (complex)
  → 逐元素乘 freqs_cis [1, L, 1, head_dim//2]
  → view_as_real → flatten → [B, L, n_heads, head_dim]
```

RoPE 只作用于 Q 和 K，不作用于 V，与原论文一致。`head_dim` 必须为偶数（有 assert 保护）。

### 9.4 FiLM Conditioning

每个 `_EncoderLayer` 包含两个 `_FiLM` 模块，分别作用于 attention 子层和 FFN 子层的 LayerNorm 输出：

```
γ, β = Linear(species_emb, 2 × d_model)   # [B, d_model] each
x_out = (1 + γ).unsqueeze(1) * LayerNorm(x) + β.unsqueeze(1)
```

`1 + γ` 的残差初始化保证训练初期 FiLM 接近 identity（`Linear` 权重和 bias 零初始化，γ ≈ 0，β ≈ 0），不破坏训练稳定性。

相比输入层 concat，FiLM 的优势：
- conv tower 直接接 4 通道 one-hot，不压缩特征维度
- 物种信息在每一层主动调制特征，而非靠 attention 被动传播
- 参数量增加：`n_layers × 2 × Linear(species_emb_dim, 2 × d_model)`（默认约 40 万）

### 9.5 AttentionPooling

用可学习的线性查询对序列维度做加权求和，将变长序列压缩为固定维度向量：

```
w = softmax(Linear(x, 1), dim=L)   # [B, L, 1]
pooled = sum(w * x, dim=L)          # [B, d_model]
```

相比 CLS token 或平均池化，对序列中信息密度不均匀的情况更鲁棒。

### 9.6 预测头（_SharedHead）

两个任务共享第一层隐藏层，然后分叉：

```
pooled [B, d_model]
    │
Linear(d_model, hidden) + GELU + Dropout   ← 共享层
    │
    ├─ phase_out: Linear(hidden, 3) → Softplus → phase_pred [B, 3]
    │                                             (ES/MS/LS log1p TPM，非负)
    │
    └─ class_logits 由 phase_pred 推算：
         WRT = (0.5×MS + LS) / (ES+MS+LS+ε)
         E_logit = (1/3 - WRT) × 10
         M_logit = -|WRT - 0.5| × 10
         L_logit = (WRT - 2/3) × 10
```

共享层强迫两个任务共享表示，有正则化效果。class 从 phase 推算保证两者预测一致，不会出现 phase 预测为 ES 主导但 class 预测为 L 的矛盾。

### 9.7 MultiTaskLoss

```
total = λ_phase × SmoothL1(phase_pred, phase_labels)
      + λ_class × CrossEntropy(class_logits, rt_class)
```

默认 `λ_phase=1.0`，`λ_class=0.5`。`class_weights` 可传入以处理类别不平衡。早停以 `val_loss_total` 为准（越低越好）。
