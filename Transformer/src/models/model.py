import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


# ── Conv block (BN → GELU → Conv, optional pool) ─────────────────────────────
class _ConvBlock(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, kernel_size: int,
                 dilation: int = 1, pool_size: int = 1,
                 dropout: float = 0.0, bn_momentum: float = 0.1,
                 norm_gamma_zero: bool = False):
        super().__init__()
        self.bn = nn.BatchNorm1d(in_ch, momentum=bn_momentum)
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
        x = self.conv(self.act(self.bn(x)))
        x = self.drop(x)
        x = self.pool(x)
        return x


# ── Dilated residual block ────────────────────────────────────────────────────
class _DilatedResidual(nn.Module):
    """Bottleneck residual block: in_ch → filters → in_ch, matching original Basenji2."""

    def __init__(self, in_ch: int, filters: int, kernel_size: int,
                 dilation: int, dropout: float, bn_momentum: float):
        super().__init__()
        self.conv1 = _ConvBlock(in_ch, filters, kernel_size,
                                dilation=dilation, bn_momentum=bn_momentum)
        # second conv: projects back to in_ch, gamma init zeros (InitZero trick)
        self.conv2 = _ConvBlock(filters, in_ch, 1,
                                dropout=dropout, bn_momentum=bn_momentum,
                                norm_gamma_zero=True)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.conv2(self.conv1(x))


# ── Basenji2 trunk ────────────────────────────────────────────────────────────
class _Basenji2Trunk(nn.Module):
    """
    Mirrors params_rt.json trunk exactly, translated to PyTorch:
      1. conv_block:       288 filters, k=15, pool=2
      2. conv_tower:       6 layers, filters 339→768 (×1.1776), k=5, pool=2
      3. dilated_residual: 11 layers, bottleneck 768→384→768, rate_mult=1.5, dropout=0.3
      4. Cropping1D:       crop 16 on each side
      5. conv_block:       1536 filters, dropout=0.05

    Input:  [B, 4, L]      (one-hot, L=32768)
    Output: [B, 1536, 224] (per-bin features)
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
                dropout=0.3, bn_momentum=bn_momentum,
            ))
            rate *= 1.5
        self.dilated = nn.Sequential(*dil_blocks)

        # 4. cropping (16 on each side) — handled in forward
        self.crop = 16

        # 5. bottleneck conv (768 → 1536)
        self.bottleneck = _ConvBlock(768, 1536, kernel_size=1, dropout=0.05, bn_momentum=bn_momentum)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: [B, 4, L]
        x = self.stem(x)
        x = self.tower(x)
        x = self.dilated(x)
        # crop 16 bins on each side → [B, 768, 224]
        x = x[:, :, self.crop:-self.crop]
        # bottleneck → [B, 1536, 224]
        return self.bottleneck(x)


# ── Prediction Heads ──────────────────────────────────────────────────────────
class _SharedHead(nn.Module):
    def __init__(self, d_model: int, hidden: int, dropout: float):
        super().__init__()
        self.shared = nn.Sequential(nn.Linear(d_model, hidden), nn.GELU(), nn.Dropout(dropout))
        self.phase_out = nn.Linear(hidden, 4)  # ES, MS, LS, G1

    def forward(self, x: torch.Tensor):
        h = self.shared(x)
        phase = F.softplus(self.phase_out(h))
        return phase


# ── Basenji2Model ─────────────────────────────────────────────────────────────
class Basenji2Model(nn.Module):
    # center 8 bins in cropped space → 1024bp centered on the 32kb window
    # 32768bp / 128bp/bin = 256 bins; crop 16 each side → 224 bins
    # center bin in cropped space = 256//2 - 16 = 112; ±4 bins → [108, 116)
    _CENTER_START = 108
    _CENTER_END   = 116

    def __init__(
        self,
        bn_momentum: float = 0.1,
        head_hidden_dim: int = 128,
        head_dropout: float = 0.1,
    ):
        super().__init__()
        self.trunk = _Basenji2Trunk(bn_momentum=bn_momentum)
        self.head = _SharedHead(1536, head_hidden_dim, head_dropout)

    def forward(self, one_hot: torch.Tensor, per_position: bool = False):
        # one_hot: [B, 4, L]
        x = self.trunk(one_hot)   # [B, 1536, 224]

        if per_position:
            # Inference: run head on every bin independently → [B, 224, 4]
            B, C, T = x.shape
            x_flat = x.permute(0, 2, 1).reshape(B * T, C)
            phase_pred = self.head(x_flat).reshape(B, T, -1)
        else:
            # Training: center 8 bins → mean → [B, 1536] → head
            x_center = x[:, :, self._CENTER_START:self._CENTER_END]  # [B, 1536, 8]
            x_pooled = x_center.mean(dim=2)                           # [B, 1536]
            phase_pred = self.head(x_pooled)                          # [B, 3]

        return {"phase_pred": phase_pred}


# ── Loss ──────────────────────────────────────────────────────────────────────
class PhaseLoss(nn.Module):
    """Poisson NLL on TPM-scale targets: loss = pred - true * log(pred + eps)."""
    def __init__(self, eps: float = 1e-6):
        super().__init__()
        self.eps = eps

    def forward(self, outputs: dict, batch: dict) -> dict:
        pred = outputs["phase_pred"]
        true = batch["phase_labels"]
        loss = (pred - true * torch.log(pred + self.eps)).mean()
        return {"total": loss, "phase": loss}
