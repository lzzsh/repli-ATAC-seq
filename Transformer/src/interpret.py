import random
import numpy as np
import torch
from typing import Callable

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
BASES = list("ACGT")


def _wrt(phase_pred: np.ndarray, eps: float = 1e-6) -> float:
    tpm = np.expm1(phase_pred)
    return float((0.5 * tpm[1] + tpm[2]) / (tpm.sum() + eps))


def _scramble(seq: str, seed: int = 42) -> str:
    rng = random.Random(seed)
    b = list(seq)
    rng.shuffle(b)
    return "".join(b)


def _gc_matched_random(seq: str, seed: int = 42) -> str:
    rng = random.Random(seed)
    gc = sum(1 for b in seq if b in "GC")
    pool = ["G", "C"] * gc + ["A", "T"] * (len(seq) - gc)
    rng.shuffle(pool)
    return "".join(pool)


# ── saturation mutagenesis ────────────────────────────────────────────────────
def saturation_mutagenesis(
    seq: str,
    predict_fn: Callable[[str], np.ndarray],
    ref_pred: np.ndarray | None = None,
) -> list[dict]:
    """Mutate every position to all 3 alt bases; return delta_ES/MS/LS/WRT per mutation."""
    if ref_pred is None:
        ref_pred = predict_fn(seq)
    ref_wrt = _wrt(ref_pred)
    results = []
    for pos, ref_base in enumerate(seq):
        for alt in BASES:
            if alt == ref_base:
                continue
            alt_pred = predict_fn(seq[:pos] + alt + seq[pos + 1:])
            results.append({
                "position": pos, "ref_base": ref_base, "alt_base": alt,
                "delta_ES": float(alt_pred[0] - ref_pred[0]),
                "delta_MS": float(alt_pred[1] - ref_pred[1]),
                "delta_LS": float(alt_pred[2] - ref_pred[2]),
                "delta_WRT": float(_wrt(alt_pred) - ref_wrt),
            })
    return results


def importance_per_position(results: list[dict]) -> np.ndarray:
    pos_max: dict[int, float] = {}
    for r in results:
        pos_max[r["position"]] = max(pos_max.get(r["position"], 0.0), abs(r["delta_WRT"]))
    arr = np.zeros(max(pos_max) + 1)
    for p, v in pos_max.items():
        arr[p] = v
    return arr


# ── motif mutation ────────────────────────────────────────────────────────────
def motif_mutation(
    seq: str, motif_start: int, motif_end: int,
    predict_fn: Callable[[str], np.ndarray],
    mode: str = "scramble",
) -> dict:
    """mode: 'scramble' | 'gc_matched_random'"""
    ref_pred = predict_fn(seq)
    motif = seq[motif_start:motif_end]
    rep = _scramble(motif) if mode == "scramble" else _gc_matched_random(motif)
    alt_pred = predict_fn(seq[:motif_start] + rep + seq[motif_end:])
    return {
        "motif_start": motif_start, "motif_end": motif_end, "mode": mode,
        "delta_ES": float(alt_pred[0] - ref_pred[0]),
        "delta_MS": float(alt_pred[1] - ref_pred[1]),
        "delta_LS": float(alt_pred[2] - ref_pred[2]),
        "delta_WRT": float(_wrt(alt_pred) - _wrt(ref_pred)),
    }


# ── flank control ─────────────────────────────────────────────────────────────
def flank_control(
    seq: str, motif_start: int, motif_end: int,
    predict_fn: Callable[[str], np.ndarray],
    flank_sizes: tuple[int, ...] = (50, 200),
) -> list[dict]:
    ref_wrt = _wrt(predict_fn(seq))
    results = []

    # motif only
    rep = _scramble(seq[motif_start:motif_end])
    p = predict_fn(seq[:motif_start] + rep + seq[motif_end:])
    results.append({"mutation_type": "motif_only", "delta_WRT": float(_wrt(p) - ref_wrt)})

    for flank in flank_sizes:
        ls, re = max(0, motif_start - flank), min(len(seq), motif_end + flank)

        # flank only (preserve motif)
        left = seq[:ls] + _scramble(seq[ls:motif_start]) + seq[motif_start:]
        both = left[:motif_end] + _scramble(left[motif_end:re]) + left[re:]
        p = predict_fn(both)
        results.append({"mutation_type": f"flank_{flank}_only", "delta_WRT": float(_wrt(p) - ref_wrt)})

        # motif + flank double
        double = seq[:ls] + _scramble(seq[ls:re]) + seq[re:]
        p = predict_fn(double)
        results.append({"mutation_type": f"motif_plus_flank_{flank}", "delta_WRT": float(_wrt(p) - ref_wrt)})

    return results


# ── Integrated Gradients ──────────────────────────────────────────────────────
def integrated_gradients(
    model,
    input_ids: torch.Tensor,
    species_id: torch.Tensor,
    target_idx: int = 3,   # 0=ES, 1=MS, 2=LS, 3=WRT
    n_steps: int = 50,
    device: str = "cpu",
) -> np.ndarray:
    """Attribution per token position [L] from [UNK] baseline."""
    model.eval()
    baseline = torch.ones_like(input_ids)
    ref_emb = model.token_emb(baseline.to(device)).detach()
    inp_emb = model.token_emb(input_ids.to(device)).detach()
    total_grads = torch.zeros_like(inp_emb)

    for step in range(1, n_steps + 1):
        interp = (ref_emb + (step / n_steps) * (inp_emb - ref_emb)).requires_grad_(True)
        sp = species_id.to(device)
        sp_emb = model.species_emb(sp).unsqueeze(1).expand(interp.shape[0], interp.shape[1], -1)
        x = model.input_proj(torch.cat([interp, sp_emb], dim=-1))
        for layer in model.layers:
            x = layer(x)
        x = model.norm(x)
        pooled = model.pooling(x)
        phase = model.phase_head(pooled)
        if target_idx < 3:
            score = phase[:, target_idx].sum()
        else:
            tpm = torch.expm1(phase)
            score = (0.5 * tpm[:, 1] + tpm[:, 2]).sum() / (tpm.sum(dim=1) + 1e-6).sum()
        score.backward()
        total_grads += interp.grad.detach()

    ig = (inp_emb - ref_emb) * total_grads / n_steps
    return ig.norm(dim=-1).squeeze(0).cpu().numpy()
