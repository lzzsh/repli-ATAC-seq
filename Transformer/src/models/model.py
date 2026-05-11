import math
import torch
import torch.nn as nn
import torch.nn.functional as F


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
    return probs / probs.amax(dim=0, keepdim=True)

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
    Enformer conv trunk (exact channel schedule from DeepMind):
      stem:   Conv1d(4→768, k=15) + Residual(ConvBlock(768)) + AttnPool/2
      tower:  6× [ConvBlock(dim_in→dim_out, k=5) + Residual(ConvBlock(dim_out)) + AttnPool/2]
              filters: 768→768→896→1024→1152→1280→1536

    Input:  [B, 4, 196608]
    Output: [B, 1536, 1536]
    """

    # exponential_linspace_int(768, 1536, num=6, divisible_by=128) prepended with 768
    _FILTER_LIST = [768, 768, 896, 1024, 1152, 1280, 1536]

    def __init__(self, bn_momentum: float = 0.1):
        super().__init__()

        # stem: plain Conv1d (no BN on raw input), then Residual(ConvBlock) + AttnPool
        self.stem_conv = nn.Conv1d(4, 768, kernel_size=15, padding=7)
        self.stem_res  = _ConvBlock(768, 768, kernel_size=1, bn_momentum=bn_momentum)
        self.stem_pool = _AttnPool(768, pool_size=2)

        # conv tower
        tower_convs, tower_res, tower_pools = [], [], []
        for dim_in, dim_out in zip(self._FILTER_LIST[:-1], self._FILTER_LIST[1:]):
            tower_convs.append(_ConvBlock(dim_in, dim_out, kernel_size=5, bn_momentum=bn_momentum))
            tower_res.append(_ConvBlock(dim_out, dim_out, kernel_size=1, bn_momentum=bn_momentum))
            tower_pools.append(_AttnPool(dim_out, pool_size=2))
        self.tower_convs = nn.ModuleList(tower_convs)
        self.tower_res   = nn.ModuleList(tower_res)
        self.tower_pools = nn.ModuleList(tower_pools)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = self.stem_conv(x)
        x = self.stem_pool(y + self.stem_res(y))
        for conv, res, pool in zip(self.tower_convs, self.tower_res, self.tower_pools):
            y = conv(x)
            x = pool(y + res(y))
        return x                                    # [B, 1536, 1536]


# ── Transformer tower ─────────────────────────────────────────────────────────
class _TransformerTower(nn.Module):
    """
    Enformer Transformer tower: n_layers× attention over full 1536 tokens.

    Input:  [B, 1536, 1536]   (C, T)
    Output: [B, 1536, 1536]   (T, C)
    """
    def __init__(self, d_model: int = 1536, n_heads: int = 8, dim_key: int = 64,
                 n_layers: int = 11, dropout: float = 0.4,
                 attn_dropout: float = 0.05, pos_dropout: float = 0.01):
        super().__init__()
        self.layers = nn.ModuleList([
            _TransformerBlock(d_model, n_heads, dim_key, dropout, attn_dropout, pos_dropout)
            for _ in range(n_layers)
        ])

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.permute(0, 2, 1)     # [B, C, T] → [B, T, C]
        for layer in self.layers:
            x = layer(x)
        return x                    # [B, T, C]


# ── EnformerModel ─────────────────────────────────────────────────────────────
class Basenji2Model(nn.Module):
    """
    Full Enformer-aligned model:
      trunk (conv)      → [B, 1536, 1536]
      transformer       → [B, 1536, 1536]  (over full 1536 tokens)
      crop 320          → [B, 896, 1536]
      final_pointwise   BN→GELU→Conv1d(1536→3072)→Dropout(0.05)→GELU → [B, 896, 3072]
      head Linear(3072→4)                                              → [B, 896, 4]
    """
    _CROP = 320

    def __init__(self, bn_momentum: float = 0.1):
        super().__init__()
        self.trunk = _EnformerTrunk(bn_momentum=bn_momentum)
        self.transformer = _TransformerTower()
        self.final_pointwise = nn.Sequential(
            nn.BatchNorm1d(1536, momentum=bn_momentum),
            nn.Conv1d(1536, 3072, kernel_size=1),
            nn.Dropout(0.05),
        )
        self.head = nn.Linear(3072, 4)

    def forward(self, one_hot: torch.Tensor):
        x = self.trunk(one_hot)                             # [B, 1536, 1536]
        x = self.transformer(x)                             # [B, 1536, 1536] (T, C)
        x = x[:, self._CROP:-self._CROP, :]                 # [B, 896, 1536]
        x = x.permute(0, 2, 1)                              # [B, 1536, 896] for BN/Conv1d
        x = _enformer_gelu(self.final_pointwise(x))         # [B, 3072, 896]
        x = x.permute(0, 2, 1)                              # [B, 896, 3072]
        return {"rt_logits": self.head(x)}                  # [B, 896, 4]


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
