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


# ── Encoder Layer (pre-norm) ──────────────────────────────────────────────────
class _EncoderLayer(nn.Module):
    def __init__(self, d_model: int, n_heads: int, dim_ff: int, dropout: float, attn_dropout: float):
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.attn = nn.MultiheadAttention(d_model, n_heads, dropout=attn_dropout, batch_first=True)
        self.ff = nn.Sequential(
            nn.Linear(d_model, dim_ff), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(dim_ff, d_model), nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.norm1(x)
        x = x + self.attn(h, h, h)[0]
        return x + self.ff(self.norm2(x))


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
        tok_dim = d_model - species_emb_dim
        self.token_emb = nn.Embedding(vocab_size, tok_dim, padding_idx=0)
        self.species_emb = nn.Embedding(n_species, species_emb_dim)
        self.input_proj = nn.Linear(d_model, d_model)
        self.layers = nn.ModuleList([
            _EncoderLayer(d_model, n_heads, dim_feedforward, dropout, attn_dropout)
            for _ in range(n_layers)
        ])
        self.norm = nn.LayerNorm(d_model)
        self.pooling = AttentionPooling(d_model)
        self.phase_head = _head(d_model, 3, head_hidden_dim, head_dropout)   # ES/MS/LS
        self.class_head = _head(d_model, 3, head_hidden_dim, head_dropout)   # E/M/L logits
        self.register_buffer(
            "freqs_cis",
            _precompute_freqs_cis(tok_dim // n_heads, max_seq_len),
            persistent=False,
        )

    def forward(self, input_ids: torch.Tensor, species_id: torch.Tensor, return_attn_weights: bool = False):
        B, L = input_ids.shape
        tok = self.token_emb(input_ids)
        sp = self.species_emb(species_id).unsqueeze(1).expand(B, L, -1)
        x = self.input_proj(torch.cat([tok, sp], dim=-1))
        for layer in self.layers:
            x = layer(x)
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
