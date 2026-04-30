import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


# ── RoPE ─────────────────────────────────────────────────────────────────────
def _precompute_freqs_cis(dim: int, seq_len: int, theta: float = 10000.0) -> torch.Tensor:
    freqs = 1.0 / (theta ** (torch.arange(0, dim, 2).float() / dim))
    t = torch.arange(seq_len, device=freqs.device)
    freqs = torch.outer(t, freqs)
    return torch.polar(torch.ones_like(freqs), freqs)  # complex64


def _apply_rope(x: torch.Tensor, freqs_cis: torch.Tensor) -> torch.Tensor:
    # x: [B, L, n_heads, head_dim]  freqs_cis: [L, head_dim//2]  (complex)
    x_ = torch.view_as_complex(x.float().reshape(*x.shape[:-1], -1, 2))
    x_rot = torch.view_as_real(x_ * freqs_cis.unsqueeze(0).unsqueeze(2)).flatten(3)
    return x_rot.type_as(x)


# ── Attention Pooling ─────────────────────────────────────────────────────────
class AttentionPooling(nn.Module):
    def __init__(self, d_model: int):
        super().__init__()
        self.query = nn.Linear(d_model, 1, bias=False)

    def forward(self, x: torch.Tensor, return_weights: bool = False):
        w = F.softmax(self.query(x), dim=1)   # [B, L, 1]
        pooled = (w * x).sum(dim=1)            # [B, D]
        return (pooled, w.squeeze(-1)) if return_weights else pooled


# ── Prediction Heads ──────────────────────────────────────────────────────────
class _SharedHead(nn.Module):
    def __init__(self, d_model: int, hidden: int, dropout: float):
        super().__init__()
        self.shared = nn.Sequential(nn.Linear(d_model, hidden), nn.GELU(), nn.Dropout(dropout))
        self.phase_out = nn.Linear(hidden, 3)
        # class logits derived from phase via WRT: no separate class_out

    def forward(self, x: torch.Tensor):
        h = self.shared(x)
        phase = F.softplus(self.phase_out(h))   # non-negative log1p TPM
        # WRT = (0.5*MS + LS) / (ES + MS + LS + eps), then threshold → E/M/L logits
        eps = 1e-6
        es, ms, ls = phase[:, 0], phase[:, 1], phase[:, 2]
        wrt = (0.5 * ms + ls) / (es + ms + ls + eps)   # [B]
        # soft logits: distance from thresholds 1/3 and 2/3
        e_logit = (1/3 - wrt) * 10
        l_logit = (wrt - 2/3) * 10
        m_logit = -torch.abs(wrt - 0.5) * 10
        class_logits = torch.stack([e_logit, m_logit, l_logit], dim=-1)  # [B, 3]
        return phase, class_logits


# ── Conv Tower (Basenji2-style) ───────────────────────────────────────────────
class _ConvStem(nn.Module):
    """Three-layer conv tower with progressive downsampling (8x total).

    Layer 1: kernel=15, stride=2, filters=128  → L/2
    Layer 2: kernel=9,  stride=2, filters=192  → L/4
    Layer 3: kernel=7,  stride=2, filters=d_model → L/8

    Input:  [B, L, d_model]   (e.g. L=4095 with stride-2 tokenizer)
    Output: [B, L//8, d_model] (e.g. ~511)
    """
    def __init__(self, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.layers = nn.ModuleList([
            self._block(d_model, 128,     kernel=15, stride=2, dropout=dropout),
            self._block(128,     192,     kernel=9,  stride=2, dropout=dropout),
            self._block(192,     d_model, kernel=7,  stride=2, dropout=dropout),
        ])

    @staticmethod
    def _block(in_ch, out_ch, kernel, stride, dropout):
        padding = kernel // 2
        return nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=kernel, stride=stride, padding=padding),
            nn.LayerNorm(out_ch),  # applied after transpose inside forward
            nn.GELU(),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: [B, L, d_model]
        x = x.transpose(1, 2)   # [B, C, L]
        for block in self.layers:
            conv, norm, act, drop = block[0], block[1], block[2], block[3]
            x = conv(x)                          # [B, C_out, L']
            x = act(norm(x.transpose(1, 2)))     # LayerNorm on last dim
            x = drop(x).transpose(1, 2)          # back to [B, C_out, L']
        return x.transpose(1, 2)                 # [B, L//8, d_model]


# ── FiLM Conditioning ─────────────────────────────────────────────────────────
class _FiLM(nn.Module):
    def __init__(self, cond_dim: int, d_model: int):
        super().__init__()
        self.proj = nn.Linear(cond_dim, 2 * d_model)
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)

    def forward(self, x: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        # x: [B, L, d_model]  cond: [B, cond_dim]
        gamma, beta = self.proj(cond).chunk(2, dim=-1)  # each [B, d_model]
        return (1 + gamma.unsqueeze(1)) * x + beta.unsqueeze(1)


# ── Encoder Layer (pre-norm) ──────────────────────────────────────────────────
class _EncoderLayer(nn.Module):
    def __init__(self, d_model: int, n_heads: int, dim_ff: int, dropout: float, attn_dropout: float, species_emb_dim: int):
        super().__init__()
        assert d_model % n_heads == 0
        self.n_heads = n_heads
        self.head_dim = d_model // n_heads
        assert self.head_dim % 2 == 0, "head_dim must be even for RoPE"
        self.scale = self.head_dim ** -0.5
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.film1 = _FiLM(species_emb_dim, d_model)
        self.film2 = _FiLM(species_emb_dim, d_model)
        self.q_proj = nn.Linear(d_model, d_model, bias=False)
        self.k_proj = nn.Linear(d_model, d_model, bias=False)
        self.v_proj = nn.Linear(d_model, d_model, bias=False)
        self.out_proj = nn.Linear(d_model, d_model, bias=False)
        self.attn_drop = nn.Dropout(attn_dropout)
        self.ff = nn.Sequential(
            nn.Linear(d_model, dim_ff), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(dim_ff, d_model), nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor, freqs_cis: torch.Tensor, species_emb: torch.Tensor, key_padding_mask: torch.Tensor | None = None) -> torch.Tensor:
        B, L, D = x.shape
        h = self.film1(self.norm1(x), species_emb)
        q = self.q_proj(h).view(B, L, self.n_heads, self.head_dim)
        k = self.k_proj(h).view(B, L, self.n_heads, self.head_dim)
        v = self.v_proj(h).view(B, L, self.n_heads, self.head_dim)
        q = _apply_rope(q, freqs_cis)
        k = _apply_rope(k, freqs_cis)
        q = q.transpose(1, 2)
        k = k.transpose(1, 2)
        v = v.transpose(1, 2)
        scores = (q @ k.transpose(-2, -1)) * self.scale
        if key_padding_mask is not None:
            scores = scores.masked_fill(key_padding_mask[:, None, None, :], float("-inf"))
        attn = torch.softmax(scores, dim=-1)
        attn = self.attn_drop(attn)
        out = (attn @ v).transpose(1, 2).reshape(B, L, D)
        x = x + self.out_proj(out)
        return x + self.ff(self.film2(self.norm2(x), species_emb))


# ── DNATransformer ────────────────────────────────────────────────────────────
class DNATransformer(nn.Module):
    def __init__(
        self,
        vocab_size: int,
        n_species: int,
        d_model: int = 256,
        n_layers: int = 6,
        n_heads: int = 8,
        dim_feedforward: int = 1024,
        dropout: float = 0.1,
        attn_dropout: float = 0.1,
        species_emb_dim: int = 64,
        head_hidden_dim: int = 128,
        head_dropout: float = 0.1,
        max_seq_len: int = 4200,
    ):
        super().__init__()
        self.max_seq_len = max_seq_len
        self.token_emb = nn.Embedding(vocab_size, d_model, padding_idx=0)
        self.conv_stem = _ConvStem(d_model, dropout)
        self.species_emb = nn.Embedding(n_species, species_emb_dim)
        self.layers = nn.ModuleList([
            _EncoderLayer(d_model, n_heads, dim_feedforward, dropout, attn_dropout, species_emb_dim)
            for _ in range(n_layers)
        ])
        self.norm = nn.LayerNorm(d_model)
        self.pooling = AttentionPooling(d_model)
        self.head = _SharedHead(d_model, head_hidden_dim, head_dropout)
        self.register_buffer(
            "freqs_cis",
            _precompute_freqs_cis(d_model // n_heads, max_seq_len // 8 + 1),
            persistent=False,
        )

    def forward(self, input_ids: torch.Tensor, species_id: torch.Tensor, key_padding_mask: torch.Tensor | None = None, return_attn_weights: bool = False):
        assert input_ids.shape[1] <= self.max_seq_len, \
            f"Input length {input_ids.shape[1]} exceeds max_seq_len {self.max_seq_len}"
        B, L = input_ids.shape
        x = self.token_emb(input_ids)                 # [B, L, d_model]
        x = self.conv_stem(x)                         # [B, L//2, d_model]
        L2 = x.shape[1]
        sp = self.species_emb(species_id)             # [B, species_emb_dim]
        freqs_cis = self.freqs_cis[:L2]
        if key_padding_mask is not None:
            key_padding_mask = key_padding_mask[:, ::8][:, :L2]  # match conv tower 8x downsample
        for layer in self.layers:
            x = layer(x, freqs_cis, sp, key_padding_mask)
        x = self.norm(x)
        if return_attn_weights:
            pooled, attn_w = self.pooling(x, return_weights=True)
        else:
            pooled, attn_w = self.pooling(x), None
        phase_pred, class_logits = self.head(pooled)
        return {
            "phase_pred": phase_pred,
            "class_logits": class_logits,
            "attn_weights": attn_w,
        }


# ── Multi-task Loss ───────────────────────────────────────────────────────────
class MultiTaskLoss(nn.Module):
    def __init__(
        self,
        class_weights: torch.Tensor | None = None,
        lambda_phase: float = 1.0,
        lambda_class: float = 0.5,
    ):
        super().__init__()
        self.lambda_phase = lambda_phase
        self.lambda_class = lambda_class
        self.smooth_l1 = nn.SmoothL1Loss()
        if class_weights is not None:
            self.register_buffer("class_weights", class_weights)
        else:
            self.class_weights = None

    def forward(self, outputs: dict, batch: dict) -> dict:
        phase_loss = self.smooth_l1(outputs["phase_pred"], batch["phase_labels"])
        class_loss = F.cross_entropy(
            outputs["class_logits"], batch["rt_class"], weight=self.class_weights,
        )
        total = self.lambda_phase * phase_loss + self.lambda_class * class_loss
        return {"total": total, "phase": phase_loss, "class": class_loss}
