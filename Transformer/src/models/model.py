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
def _head(d_model: int, out_dim: int, hidden: int, dropout: float) -> nn.Sequential:
    return nn.Sequential(
        nn.Linear(d_model, hidden), nn.GELU(), nn.Dropout(dropout),
        nn.Linear(hidden, out_dim),
    )


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
        max_seq_len: int = 2800,
    ):
        super().__init__()
        self.token_emb = nn.Embedding(vocab_size, d_model, padding_idx=0)
        self.species_emb = nn.Embedding(n_species, species_emb_dim)
        self.layers = nn.ModuleList([
            _EncoderLayer(d_model, n_heads, dim_feedforward, dropout, attn_dropout, species_emb_dim)
            for _ in range(n_layers)
        ])
        self.norm = nn.LayerNorm(d_model)
        self.pooling = AttentionPooling(d_model)
        self.phase_head = _head(d_model, 3, head_hidden_dim, head_dropout)
        self.class_head = _head(d_model, 3, head_hidden_dim, head_dropout)
        self.register_buffer(
            "freqs_cis",
            _precompute_freqs_cis(d_model // n_heads, max_seq_len),
            persistent=False,
        )

    def forward(self, input_ids: torch.Tensor, species_id: torch.Tensor, key_padding_mask: torch.Tensor | None = None, return_attn_weights: bool = False):
        B, L = input_ids.shape
        x = self.token_emb(input_ids)                 # [B, L, d_model]
        sp = self.species_emb(species_id)             # [B, species_emb_dim]
        freqs_cis = self.freqs_cis[:L]
        for layer in self.layers:
            x = layer(x, freqs_cis, sp, key_padding_mask)
        x = self.norm(x)
        if return_attn_weights:
            pooled, attn_w = self.pooling(x, return_weights=True)
        else:
            pooled, attn_w = self.pooling(x), None
        return {
            "phase_pred": self.phase_head(pooled),
            "class_logits": self.class_head(pooled),
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
