from __future__ import annotations
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint_sequential

from ..data.dataset import SpeciesConfig
from .config_model import RepliformerConfig


# ── Positional feature functions (from Enformer) ──────────────────────────────

def _pos_features_exponential(positions, features, seq_len):
    max_range = math.log(seq_len) / math.log(2.)
    half_life = 2 ** torch.linspace(math.log(3.) / math.log(2.), max_range, features,
                                    device=positions.device)
    return torch.exp(-math.log(2.) / half_life * positions.abs().unsqueeze(-1))

def _pos_features_central_mask(positions, features, seq_len):
    center_widths = (2 ** torch.arange(1, features + 1, device=positions.device).float()) - 1
    return (center_widths.unsqueeze(0) > positions.abs().unsqueeze(-1)).float()

def _gamma_pdf(x, concentration, rate):
    log_unnorm = torch.xlogy(concentration - 1., x) - rate * x
    log_norm = torch.lgamma(concentration) - concentration * torch.log(rate)
    return torch.exp(log_unnorm - log_norm)

def _pos_features_gamma(positions, features, seq_len, eps=1e-8):
    stddev = seq_len / (2 * features)
    start_mean = seq_len / features
    mean = torch.linspace(start_mean, seq_len, features, device=positions.device)
    concentration = (mean / stddev) ** 2
    rate = mean / stddev ** 2
    probs = _gamma_pdf(positions.abs().unsqueeze(-1), concentration, rate) + eps
    return probs / probs.amax(dim=-1, keepdim=True)

def _get_positional_embed(seq_len: int, feature_size: int, device) -> torch.Tensor:
    """[2T-1, feature_size] positional basis features (Enformer-style)."""
    assert feature_size % 6 == 0, f'feature_size must be divisible by 6, got {feature_size}'
    n = feature_size // 6
    distances = torch.arange(-(seq_len - 1), seq_len, device=device, dtype=torch.float32)
    parts = [
        _pos_features_exponential(distances, n, seq_len),
        _pos_features_central_mask(distances, n, seq_len),
        _pos_features_gamma(distances, n, seq_len),
    ]
    embeddings = torch.cat(parts, dim=-1)                          # [2T-1, feature_size//2]
    sign = torch.sign(distances).unsqueeze(-1)                     # [2T-1, 1]
    return torch.cat([embeddings, sign * embeddings], dim=-1)      # [2T-1, feature_size]

def _relative_shift(x: torch.Tensor) -> torch.Tensor:
    """TransformerXL relative shift. x: [B, H, T, 2T-1] → [B, H, T, T]"""
    to_pad = torch.zeros_like(x[..., :1])
    x = torch.cat((to_pad, x), dim=-1)
    B, H, t1, t2 = x.shape
    x = x.reshape(B, H, t2, t1)
    x = x[:, :, 1:, :]
    x = x.reshape(B, H, t1, t2 - 1)
    return x[..., :((t2 + 1) // 2)]


def _enformer_gelu(x: torch.Tensor) -> torch.Tensor:
    """Enformer GELU: sigmoid(1.702*x)*x (swish approximation used throughout)."""
    return torch.sigmoid(1.702 * x) * x


# ── Conv block ────────────────────────────────────────────────────────────────
class _ConvBlock(nn.Module):
    """Enformer ConvBlock: BN → GELU → Conv1d."""
    def __init__(self, in_ch: int, out_ch: int, kernel_size: int = 1,
                 bn_momentum: float = 0.1):
        super().__init__()
        self.bn   = nn.BatchNorm1d(in_ch, momentum=bn_momentum)
        self.conv = nn.Conv1d(in_ch, out_ch, kernel_size, padding=kernel_size // 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.conv(_enformer_gelu(self.bn(x)))



# ── Attention pooling ─────────────────────────────────────────────────────────
class _AttnPool(nn.Module):
    """Per-channel softmax attention pooling (matches Enformer AttentionPool).

    Uses Conv2d(dim, dim, 1, bias=False) with dirac_ init × 2.
    Handles odd-length inputs by padding then masking.

    Input:  [B, C, L]
    Output: [B, C, L // pool_size]
    """
    def __init__(self, in_channels: int, pool_size: int = 2):
        super().__init__()
        self.pool_size = pool_size
        self.to_attn_logits = nn.Conv2d(in_channels, in_channels, 1, bias=False)
        nn.init.dirac_(self.to_attn_logits.weight)
        with torch.no_grad():
            self.to_attn_logits.weight.mul_(2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, C, L = x.shape
        remainder = L % self.pool_size
        needs_padding = remainder > 0
        if needs_padding:
            x = F.pad(x, (0, self.pool_size - remainder), value=0)
            mask = torch.zeros((B, 1, L), dtype=torch.bool, device=x.device)
            mask = F.pad(mask, (0, self.pool_size - remainder), value=True)

        # [B, C, L', p] where L' = L_padded // pool_size
        x = x.reshape(B, C, -1, self.pool_size)
        logits = self.to_attn_logits(x)

        if needs_padding:
            mask = mask.reshape(B, 1, -1, self.pool_size)
            mask_val = -torch.finfo(logits.dtype).max
            logits = logits.masked_fill(mask, mask_val)

        attn = logits.softmax(dim=-1)
        return (x * attn).sum(dim=-1)              # [B, C, L']


# ── Enformer attention ────────────────────────────────────────────────────────
class _EnformerAttention(nn.Module):
    """
    Enformer multi-head attention (matches modeling_enformer.py Attention class):
    - dim_key=64 per head (Q, K); dim_value=d_model//heads=192 per head (V)
    - TransformerXL relative position: basis features → to_rel_k
    - rel_content_bias (r_w_bias) and rel_pos_bias (r_r_bias)
    - relative_shift for causal-free relative logits
    - attn_dropout=0.05, pos_dropout=0.01
    - Q/K/V: bias=False; out_proj: zero-initialized
    """
    def __init__(self, d_model: int = 1536, n_heads: int = 8, dim_key: int = 64,
                 attn_dropout: float = 0.05, pos_dropout: float = 0.01):
        super().__init__()
        self.scale = dim_key ** -0.5
        self.n_heads = n_heads
        self.dim_key = dim_key
        dim_value = d_model // n_heads                  # 192
        num_rel_pos_features = d_model // n_heads       # 192

        self.to_q = nn.Linear(d_model, dim_key * n_heads, bias=False)
        self.to_k = nn.Linear(d_model, dim_key * n_heads, bias=False)
        self.to_v = nn.Linear(d_model, dim_value * n_heads, bias=False)

        self.to_out = nn.Linear(dim_value * n_heads, d_model)
        nn.init.zeros_(self.to_out.weight)
        nn.init.zeros_(self.to_out.bias)

        self.to_rel_k = nn.Linear(num_rel_pos_features, dim_key * n_heads, bias=False)
        self.rel_content_bias = nn.Parameter(torch.randn(1, n_heads, 1, dim_key))
        self.rel_pos_bias     = nn.Parameter(torch.randn(1, n_heads, 1, dim_key))

        self.pos_dropout  = nn.Dropout(pos_dropout)
        self.attn_dropout = nn.Dropout(attn_dropout)
        self.num_rel_pos_features = num_rel_pos_features

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, T, _ = x.shape
        H, D = self.n_heads, self.dim_key
        Dv = x.shape[-1] // H

        q = self.to_q(x).reshape(B, T, H, D).permute(0, 2, 1, 3)   # [B, H, T, D]
        k = self.to_k(x).reshape(B, T, H, D).permute(0, 2, 1, 3)
        v = self.to_v(x).reshape(B, T, H, Dv).permute(0, 2, 1, 3)  # [B, H, T, Dv]

        q = q * self.scale

        # content logits: (q + r_w_bias) @ k^T
        content_logits = torch.einsum('bhid,bhjd->bhij', q + self.rel_content_bias, k)

        # relative position logits
        positions = _get_positional_embed(T, self.num_rel_pos_features, x.device)
        positions = self.pos_dropout(positions)
        rel_k = self.to_rel_k(positions)                             # [2T-1, H*D]
        rel_k = rel_k.reshape(-1, H, D).permute(1, 0, 2)            # [H, 2T-1, D]
        rel_logits = torch.einsum('bhid,hjd->bhij', q + self.rel_pos_bias, rel_k)
        rel_logits = _relative_shift(rel_logits)                     # [B, H, T, T]

        attn = (content_logits + rel_logits).softmax(dim=-1)
        attn = self.attn_dropout(attn)

        out = torch.einsum('bhij,bhjd->bhid', attn, v)
        out = out.permute(0, 2, 1, 3).reshape(B, T, -1)             # [B, T, H*Dv]
        return self.to_out(out)


# ── Transformer block ─────────────────────────────────────────────────────────
class _TransformerBlock(nn.Module):
    """Enformer Transformer layer: pre-norm attention + pre-norm FFN with ReLU."""

    def __init__(self, d_model: int = 1536, n_heads: int = 8, dim_key: int = 64,
                 dropout: float = 0.4, attn_dropout: float = 0.05, pos_dropout: float = 0.01):
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.attn  = _EnformerAttention(d_model, n_heads, dim_key, attn_dropout, pos_dropout)
        self.drop1 = nn.Dropout(dropout)

        self.norm2 = nn.LayerNorm(d_model)
        self.ff1   = nn.Linear(d_model, d_model * 2)
        self.ff2   = nn.Linear(d_model * 2, d_model)
        self.drop2 = nn.Dropout(dropout)
        self.drop3 = nn.Dropout(dropout)
        self.act   = nn.ReLU()          # Enformer FFN uses ReLU

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.drop1(self.attn(self.norm1(x)))
        h = self.norm2(x)
        h = self.ff2(self.drop2(self.act(self.ff1(h))))
        return x + self.drop3(h)


# ── Enformer trunk ────────────────────────────────────────────────────────────
class _EnformerTrunk(nn.Module):
    """
    Enformer conv trunk.
      stem:   Conv1d(4→stem_channels, k=stem_kernel_size) + Residual(ConvBlock) + AttnPool/2
      tower:  len(tower_filter_list)-1 × [ConvBlock + Residual + AttnPool/2]
    """

    def __init__(self, cfg: RepliformerConfig):
        super().__init__()
        stem_ch = cfg.stem_channels
        filters = cfg.tower_filter_list
        bm = cfg.bn_momentum
        n_tower = len(filters) - 1

        self.stem_conv = nn.Conv1d(4, stem_ch, kernel_size=cfg.stem_kernel_size,
                                   padding=cfg.stem_kernel_size // 2)
        self.stem_res  = _ConvBlock(stem_ch, stem_ch, kernel_size=1, bn_momentum=bm)
        self.stem_pool = _AttnPool(stem_ch, pool_size=2)

        tower_convs, tower_res, tower_pools = [], [], []
        for i, (dim_in, dim_out) in enumerate(zip(filters[:-1], filters[1:])):
            tower_convs.append(_ConvBlock(dim_in, dim_out, kernel_size=5, bn_momentum=bm))
            tower_res.append(_ConvBlock(dim_out, dim_out, kernel_size=1, bn_momentum=bm))
            use_avg = (cfg.pool_type == "avg_late" and i >= n_tower - 3)
            tower_pools.append(
                nn.AvgPool1d(kernel_size=2, stride=2) if use_avg
                else _AttnPool(dim_out, pool_size=2)
            )
        self.tower_convs = nn.ModuleList(tower_convs)
        self.tower_res   = nn.ModuleList(tower_res)
        self.tower_pools = nn.ModuleList(tower_pools)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = self.stem_conv(x)
        x = self.stem_pool(y + self.stem_res(y))
        for conv, res, pool in zip(self.tower_convs, self.tower_res, self.tower_pools):
            y = conv(x)
            x = pool(y + res(y))
        return x


# ── Transformer tower ─────────────────────────────────────────────────────────
class _TransformerTower(nn.Module):
    """Enformer Transformer tower."""
    def __init__(self, cfg: RepliformerConfig):
        super().__init__()
        self.layers = nn.ModuleList([
            _TransformerBlock(cfg.d_model, cfg.n_heads, cfg.dim_key,
                              cfg.dropout, cfg.attn_dropout, cfg.pos_dropout)
            for _ in range(cfg.n_layers)
        ])

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.permute(0, 2, 1)     # [B, C, T] → [B, T, C]
        # gradient checkpointing: recompute activations during backward
        # trades compute for memory — required for T=1536 on 40GB GPUs
        x = checkpoint_sequential(self.layers, len(self.layers), x,
                                  use_reentrant=False)
        return x                    # [B, T, C]


# ── EnformerModel ─────────────────────────────────────────────────────────────
class RepliformerModel(nn.Module):
    """
    Full Enformer-aligned model driven by RepliformerConfig.
      trunk (conv)      → [B, d_model, T]
      transformer       → [B, T, d_model]
      crop cfg.crop_size → [B, n_output_bins, d_model]
      final_pointwise   BN→GELU→Conv1d(d_model→pointwise_channels)→Dropout→GELU
      _heads[head: str] Linear(pointwise_channels→n_output_bins)
    """

    def __init__(self, species_configs: list[SpeciesConfig],
                 model_cfg: RepliformerConfig | None = None):
        super().__init__()
        assert len(species_configs) > 0, "species_configs must not be empty"
        cfg = model_cfg if model_cfg is not None else RepliformerConfig()
        self._cfg = cfg
        self.trunk       = _EnformerTrunk(cfg)
        self.transformer = _TransformerTower(cfg)
        trunk_out = cfg.tower_filter_list[-1]
        self.trunk_proj = (
            nn.Linear(trunk_out, cfg.d_model, bias=False)
            if trunk_out != cfg.d_model else nn.Identity()
        )
        self.final_bn    = nn.BatchNorm1d(cfg.d_model, momentum=cfg.bn_momentum)
        self.final_conv  = nn.Conv1d(cfg.d_model, cfg.pointwise_channels, kernel_size=1)
        self.final_drop  = nn.Dropout(cfg.final_dropout)
        self._heads = nn.ModuleDict({
            sp.name: nn.Linear(cfg.pointwise_channels, cfg.n_output_bins)
            for sp in species_configs
        })

    def forward(self, one_hot: torch.Tensor, head: str) -> dict:
        crop = self._cfg.crop_size
        x = self.trunk(one_hot)                  # [B, trunk_out, T]
        x = self.trunk_proj(x.permute(0, 2, 1)).permute(0, 2, 1)  # [B, d_model, T]
        x = self.transformer(x)
        x = x[:, crop:-crop, :]
        x = x.permute(0, 2, 1)
        x = _enformer_gelu(self.final_bn(x))
        x = _enformer_gelu(self.final_drop(self.final_conv(x)))
        x = x.permute(0, 2, 1)
        return {"rt_signals": F.softplus(self._heads[head](x))}


# ── Loss ──────────────────────────────────────────────────────────────────────
class RTSignalLoss(nn.Module):
    def forward(self, outputs: dict, batch: dict) -> dict:
        pred   = outputs["rt_signals"]   # [B, T, 4]
        target = batch["rt_signals"]     # [B, T, 4]
        mask = ~torch.isnan(target)
        loss = F.smooth_l1_loss(pred[mask], target[mask])
        return {"total": loss, "rt": loss}
