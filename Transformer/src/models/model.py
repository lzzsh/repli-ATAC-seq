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
    Enformer trunk at 196608bp input:
      1. conv_block + AttnPool2: 4→288, k=15
      2. conv_tower + AttnPool2: 6 layers, filters 339→768 (×1.1776), k=5
      3. conv_block: 768→1536, k=1, dropout=0.05
      4. dilated residual tower: 11 layers, k=3, dilation=[1,2,4,8,16,32,64,128,256,512,1]
      5. Cropping1D: crop 320 on each side (at 128bp resolution)

    Input:  [B, 4, 196608]
    Output: [B, 1536, 896]  (896 bins × 128bp, after crop-320 each side)
    """

    _TOWER_FILTERS = [int(np.round(339 * (1.1776 ** i))) for i in range(6)]
    _DILATIONS     = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1]

    def __init__(self, bn_momentum: float = 0.1):
        super().__init__()

        # 1. stem: conv + attn pool
        self.stem_conv = _ConvBlock(4, 288, kernel_size=15, bn_momentum=bn_momentum)
        self.stem_pool = _AttnPool(288, pool_size=2)

        # 2. conv tower: 6 layers, each conv + attn pool
        tower_convs, tower_pools = [], []
        in_ch = 288
        for out_ch in self._TOWER_FILTERS:
            tower_convs.append(_ConvBlock(in_ch, out_ch, kernel_size=5, bn_momentum=bn_momentum))
            tower_pools.append(_AttnPool(out_ch, pool_size=2))
            in_ch = out_ch
        self.tower_convs = nn.ModuleList(tower_convs)
        self.tower_pools = nn.ModuleList(tower_pools)

        # 3. bottleneck conv (768 → 1536)
        self.bottleneck = _ConvBlock(768, 1536, kernel_size=1, dropout=0.05, bn_momentum=bn_momentum)

        # 4. dilated residual tower: 11 layers matching Enformer
        self.dilated_tower = nn.ModuleList([
            _DilatedResidual(in_ch=1536, filters=1536, kernel_size=3,
                             dilation=d, dropout=0.3, bn_momentum=bn_momentum)
            for d in self._DILATIONS
        ])

        # 5. cropping (320 on each side) — handled in forward
        self.crop = 320

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.stem_pool(self.stem_conv(x))
        for conv, pool in zip(self.tower_convs, self.tower_pools):
            x = pool(conv(x))
        x = self.bottleneck(x)              # [B, 1536, L]
        for block in self.dilated_tower:
            x = block(x)
        x = x[:, :, self.crop:-self.crop]   # crop 320 → [B, 1536, 896]
        return x


# ── Prediction Head ───────────────────────────────────────────────────────────
class _RTHead(nn.Module):
    """Transformer head: n_layers× self-attention → classify at 128bp resolution.

    Input:  [B, in_features, T]
    Output: [B, T, n_classes]
    """
    def __init__(
        self,
        in_features: int = 1536,
        d_model: int = 512,
        n_heads: int = 8,
        ffn_dim: int = 1024,
        n_layers: int = 4,
        n_classes: int = 4,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.proj = nn.Linear(in_features, d_model) if in_features != d_model else nn.Identity()
        self.pos_bias = _RelativePosBias(n_heads=n_heads, max_len=896)
        self.layers = nn.ModuleList([
            _TransformerBlock(d_model, n_heads, ffn_dim, dropout)
            for _ in range(n_layers)
        ])
        self.norm = nn.LayerNorm(d_model)
        self.out = nn.Linear(d_model, n_classes)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: [B, C, T]
        x = x.permute(0, 2, 1)             # [B, T, C]
        x = self.proj(x)                   # [B, T, d_model]
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
            in_features=1536, d_model=512, n_heads=8,
            ffn_dim=1024, n_layers=4, n_classes=4, dropout=0.2,
        )

    def forward(self, one_hot: torch.Tensor):
        x = self.trunk(one_hot)             # [B, 1536, T]
        return {"rt_logits": self.head(x)}  # [B, T, 4]


# ── Loss ──────────────────────────────────────────────────────────────────────
class RTClassLoss(nn.Module):
    def __init__(self, class_weights: list[float] | None = None):
        super().__init__()
        if class_weights is not None:
            self.register_buffer("weight", torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.weight = None

    def forward(self, outputs: dict, batch: dict) -> dict:
        logits = outputs["rt_logits"]
        labels = batch["rt_labels"]
        loss = F.cross_entropy(logits.reshape(-1, 4), labels.reshape(-1), weight=self.weight)
        return {"total": loss, "rt": loss}
