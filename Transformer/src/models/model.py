import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


# ── Conv block ────────────────────────────────────────────────────────────────
# Two modes controlled by `pre_norm`:
#   pre_norm=True  (stem/tower/bottleneck): BN → GELU → Conv  (pre-activation)
#   pre_norm=False (residual internals):    GELU → Conv → BN   (matches Basenji2 paper)
class _ConvBlock(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, kernel_size: int,
                 dilation: int = 1, pool_size: int = 1,
                 dropout: float = 0.0, bn_momentum: float = 0.1,
                 norm_gamma_zero: bool = False, pre_norm: bool = True):
        super().__init__()
        self.pre_norm = pre_norm
        self.bn = nn.BatchNorm1d(in_ch if pre_norm else out_ch, momentum=bn_momentum)
        if norm_gamma_zero:
            nn.init.zeros_(self.bn.weight)
        self.act = nn.GELU()
        self.conv = nn.Conv1d(
            in_ch, out_ch, kernel_size,
            padding=dilation * (kernel_size // 2),
            dilation=dilation, bias=False,
        )
        self.drop = nn.Dropout(dropout) if dropout > 0 else nn.Identity()
        self.pool = nn.MaxPool1d(pool_size) if pool_size > 1 else nn.Identity()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.pre_norm:
            x = self.conv(self.act(self.bn(x)))
        else:
            x = self.bn(self.conv(self.act(x)))
        x = self.drop(x)
        x = self.pool(x)
        return x


# ── Dilated residual block ────────────────────────────────────────────────────
class _DilatedResidual(nn.Module):
    """Bottleneck residual block: in_ch → filters → in_ch, matching original Basenji2."""

    def __init__(self, in_ch: int, filters: int, kernel_size: int,
                 dilation: int, dropout: float, bn_momentum: float):
        super().__init__()
        # paper order: GELU → Conv → BN (pre_norm=False)
        self.conv1 = _ConvBlock(in_ch, filters, kernel_size,
                                dilation=dilation, bn_momentum=bn_momentum,
                                pre_norm=False)
        # second conv: projects back to in_ch, gamma init zeros (InitZero trick)
        self.conv2 = _ConvBlock(filters, in_ch, 1,
                                dropout=dropout, bn_momentum=bn_momentum,
                                norm_gamma_zero=True, pre_norm=False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.conv2(self.conv1(x))


# ── Attention pooling ─────────────────────────────────────────────────────────
class _AttnPool(nn.Module):
    """Attention pooling: replaces MaxPool. Learns softmax weights per position.

    For each pool window of size `pool_size`, computes a scalar logit per position
    via a linear layer, applies softmax, then returns the weighted sum of features.

    Input:  [B, C, L]
    Output: [B, C, L // pool_size]
    """
    def __init__(self, in_channels: int, pool_size: int = 2):
        super().__init__()
        self.pool_size = pool_size
        self.weight = nn.Linear(in_channels, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, C, L = x.shape
        L_trunc = (L // self.pool_size) * self.pool_size
        x = x[:, :, :L_trunc]                              # [B, C, L_trunc]
        x = x.reshape(B, C, L_trunc // self.pool_size, self.pool_size)
        x_t = x.permute(0, 2, 3, 1)                        # [B, L//p, p, C]
        w = self.weight(x_t)                                # [B, L//p, p, 1]
        w = torch.softmax(w, dim=2)
        out = (x_t * w).sum(dim=2)                          # [B, L//p, C]
        return out.permute(0, 2, 1)                         # [B, C, L//p]


# ── Relative positional bias ──────────────────────────────────────────────────
class _RelativePosBias(nn.Module):
    """Learnable relative positional bias added to attention logits (Enformer-style).

    bias table has shape [n_heads, 2*max_len - 1].
    For positions i, j: bias[h, i, j] = table[h, j - i + max_len - 1]
    """
    def __init__(self, n_heads: int, max_len: int):
        super().__init__()
        self.max_len = max_len
        self.bias = nn.Parameter(torch.zeros(n_heads, 2 * max_len - 1))

    def forward(self, seq_len: int) -> torch.Tensor:
        idx = torch.arange(seq_len, device=self.bias.device)
        offsets = idx.unsqueeze(1) - idx.unsqueeze(0)          # [T, T]
        table_idx = offsets + self.max_len - 1                  # [T, T], 0-based
        return self.bias[:, table_idx]                          # [n_heads, T, T]



# ── Transformer block ─────────────────────────────────────────────────────────
class _TransformerBlock(nn.Module):
    """Single Transformer layer: multi-head self-attention + FFN, with relative pos bias."""

    def __init__(self, d_model: int, n_heads: int, ffn_dim: int, dropout: float = 0.1):
        super().__init__()
        self.attn = nn.MultiheadAttention(
            d_model, n_heads, dropout=dropout, batch_first=True
        )
        self.ff1 = nn.Linear(d_model, ffn_dim)
        self.ff2 = nn.Linear(ffn_dim, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.drop = nn.Dropout(dropout)
        self.act = nn.GELU()

    def forward(self, x: torch.Tensor, attn_bias: torch.Tensor) -> torch.Tensor:
        # attn_bias: [n_heads, T, T] → expand to [B*n_heads, T, T]
        B, T, _ = x.shape
        n_heads = attn_bias.shape[0]
        mask = attn_bias.unsqueeze(0).expand(B, -1, -1, -1).reshape(B * n_heads, T, T)

        # pre-norm attention
        h = self.norm1(x)
        h, _ = self.attn(h, h, h, attn_mask=mask, need_weights=False)
        x = x + self.drop(h)

        # pre-norm FFN
        h = self.norm2(x)
        h = self.ff2(self.drop(self.act(self.ff1(h))))
        return x + self.drop(h)


# ── Enformer trunk ────────────────────────────────────────────────────────────
class _EnformerTrunk(nn.Module):
    """
    Enformer-style trunk at 196608bp input:
      1. conv_block:       288 filters, k=15, pool=2
      2. conv_tower:       6 layers, filters 339→768 (×1.1776), k=5, pool=2
      3. dilated_residual: 11 layers, bottleneck 768→384→768, rate_mult=1.5, dropout=0.4
      4. Cropping1D:       crop 16 on each side  (at 128bp resolution)
      5. conv_block:       3072 filters, dropout=0.05

    Input:  [B, 4, 196608]
    Output: [B, 3072, 1504]  (1504 bins × 128bp, after crop-16 each side)
    """

    # conv tower filter schedule: 339 * 1.1776^i, rounded, 6 layers
    _TOWER_FILTERS = [int(np.round(339 * (1.1776 ** i))) for i in range(6)]  # [339,399,470,554,652,768]

    def __init__(self, bn_momentum: float = 0.1):
        super().__init__()

        # 1. stem conv
        self.stem = _ConvBlock(4, 288, kernel_size=15, pool_size=2, bn_momentum=bn_momentum)

        # 2. conv tower (6 layers, each pool=2 → total 64x downsample with stem)
        tower = []
        in_ch = 288
        for out_ch in self._TOWER_FILTERS:
            tower.append(_ConvBlock(in_ch, out_ch, kernel_size=5, pool_size=2, bn_momentum=bn_momentum))
            in_ch = out_ch
        self.tower = nn.Sequential(*tower)

        # 3. dilated residual tower (11 layers, bottleneck 768→384→768, rate_mult=1.5)
        dil_blocks = []
        rate = 1.0
        for _ in range(11):
            dil_blocks.append(_DilatedResidual(
                in_ch=768, filters=384, kernel_size=3,
                dilation=int(np.round(rate)),
                dropout=0.4, bn_momentum=bn_momentum,
            ))
            rate *= 1.5
        self.dilated = nn.Sequential(*dil_blocks)

        # 4. cropping (16 on each side) — handled in forward
        self.crop = 16

        # 5. bottleneck conv (768 → 3072)
        self.bottleneck = _ConvBlock(768, 3072, kernel_size=1, dropout=0.05, bn_momentum=bn_momentum)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: [B, 4, L]
        x = self.stem(x)
        x = self.tower(x)
        x = self.dilated(x)
        x = x[:, :, self.crop:-self.crop]  # crop 16 → [B, 768, L/128-32]
        return self.bottleneck(x)          # → [B, 3072, L/128-32]


# ── Prediction Head ───────────────────────────────────────────────────────────
class _RTHead(nn.Module):
    """Transformer head: project → 8× self-attention → classify at 128bp resolution.

    Input:  [B, in_features, T]   (T = 1504 for 196608bp input)
    Output: [B, T, n_classes]
    """
    def __init__(
        self,
        in_features: int = 3072,
        d_model: int = 512,
        n_heads: int = 8,
        ffn_dim: int = 2048,
        n_layers: int = 8,
        n_classes: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.proj = nn.Linear(in_features, d_model)
        self.pos_bias = _RelativePosBias(n_heads=n_heads, max_len=1504)
        self.layers = nn.ModuleList([
            _TransformerBlock(d_model, n_heads, ffn_dim, dropout)
            for _ in range(n_layers)
        ])
        self.norm = nn.LayerNorm(d_model)
        self.out = nn.Linear(d_model, n_classes)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: [B, C, T]
        x = x.permute(0, 2, 1)             # [B, T, C]
        x = self.proj(x)                    # [B, T, d_model]
        T = x.shape[1]
        attn_bias = self.pos_bias(T)        # [n_heads, T, T]
        for layer in self.layers:
            x = layer(x, attn_bias)
        x = self.norm(x)
        return self.out(x)                  # [B, T, n_classes]


# ── Basenji2Model ─────────────────────────────────────────────────────────────
class Basenji2Model(nn.Module):
    def __init__(self, bn_momentum: float = 0.1):
        super().__init__()
        self.trunk = _EnformerTrunk(bn_momentum=bn_momentum)
        self.head = _RTHead(
            in_features=3072, d_model=512, n_heads=8,
            ffn_dim=2048, n_layers=8, n_classes=4, dropout=0.1,
        )

    def forward(self, one_hot: torch.Tensor):
        # one_hot: [B, 4, L]
        x = self.trunk(one_hot)             # [B, 3072, T]
        return {"rt_logits": self.head(x)}  # [B, T, 4]


# ── Loss ──────────────────────────────────────────────────────────────────────
class RTClassLoss(nn.Module):
    def forward(self, outputs: dict, batch: dict) -> dict:
        logits = outputs["rt_logits"]          # [B, 28, 4]
        labels = batch["rt_labels"]            # [B, 28] long
        loss = F.cross_entropy(logits.reshape(-1, 4), labels.reshape(-1))
        return {"total": loss, "rt": loss}
