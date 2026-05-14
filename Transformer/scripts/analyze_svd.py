#!/usr/bin/env python3
"""
Layer-wise SVD activation analysis for RepliformerModel.

Usage:
    python scripts/analyze_svd.py \
        --checkpoint /liaozizhuo/repli-ATAC-seq/outputs/basenji2_wt/checkpoints/best_model.pt \
        --config src/configs/transformer_wt.yaml \
        --species rice \
        --n_samples 32 \
        --top_k 64 \
        --out_dir outputs/svd_analysis
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import yaml
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.data.dataset import load_manifest, RepliSeqDataset
from src.models.model import RepliformerModel
from src.models.config_model import RepliformerConfig


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True, help="Path to .pt checkpoint")
    p.add_argument("--config",     required=True, help="Path to training YAML config")
    p.add_argument("--species",    default="rice", help="Species head to use")
    p.add_argument("--n_samples",  type=int, default=32, help="Samples for forward pass")
    p.add_argument("--top_k",      type=int, default=64,  help="Singular values to plot")
    p.add_argument("--out_dir",    default="outputs/svd_analysis")
    return p.parse_args()


def load_model(checkpoint_path: str, config_path: str, device: torch.device):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    species_configs = load_manifest(cfg["data"]["manifest"])
    model = RepliformerModel(
        species_configs=species_configs,
        model_cfg=RepliformerConfig(**cfg.get("model", {})),
    ).to(device)
    ckpt = torch.load(checkpoint_path, map_location=device)
    state = ckpt["model"] if "model" in ckpt else ckpt
    model.load_state_dict(state)
    model.eval()
    return model, species_configs, cfg


def load_data(species_configs, species_name: str, cfg: dict,
              n_samples: int, device: torch.device) -> torch.Tensor:
    sp = next(s for s in species_configs if s.name == species_name)
    ds = RepliSeqDataset(
        [sp], split="val",
        window_size=cfg["data"]["input_window_length"],
        rc_prob=0.0, shift_max=0,
    )
    loader = DataLoader(ds, batch_size=1, shuffle=False, num_workers=0)
    one_hots = []
    for batch in loader:
        one_hots.append(batch["one_hot"])
        if len(one_hots) >= n_samples:
            break
    return torch.cat(one_hots, dim=0).to(device)   # [N, 4, L]


def register_hooks(model: RepliformerModel) -> Tuple[Dict[str, List], List]:
    activations: Dict[str, List] = defaultdict(list)
    handles = []

    def make_hook(name: str):
        def hook(module, input, output):
            if isinstance(output, torch.Tensor):
                activations[name].append(output.detach().cpu())
        return hook

    trunk = model.trunk

    handles.append(trunk.stem_conv.register_forward_hook(make_hook("stem_conv")))
    handles.append(trunk.stem_res.register_forward_hook(make_hook("stem_res")))
    handles.append(trunk.stem_pool.register_forward_hook(make_hook("stem_pool")))

    for i, (conv, res, pool) in enumerate(
        zip(trunk.tower_convs, trunk.tower_res, trunk.tower_pools)
    ):
        handles.append(conv.register_forward_hook(make_hook(f"tower_conv_{i}")))
        handles.append(res.register_forward_hook(make_hook(f"tower_res_{i}")))
        handles.append(pool.register_forward_hook(make_hook(f"tower_pool_{i}")))

    handles.append(model.trunk_proj.register_forward_hook(make_hook("trunk_proj")))

    for i, block in enumerate(model.transformer.layers):
        handles.append(block.attn.register_forward_hook(make_hook(f"transformer_{i}_attn_out")))
        handles.append(block.register_forward_hook(make_hook(f"transformer_{i}_block_out")))

    handles.append(model.final_bn.register_forward_hook(make_hook("final_bn")))
    handles.append(model.final_conv.register_forward_hook(make_hook("final_conv")))

    return activations, handles


def _reshape_activation(name: str, act: torch.Tensor) -> torch.Tensor:
    """
    Conv layers  [B, C, L] → [B*L, C]  (feature axis = C)
    Transformer  [B, T, D] → [B*T, D]  (feature axis = D)
    Both cases: act.ndim==3, feature dim is dim-1 for transformer, dim-2 for conv.
    We distinguish by name prefix.
    """
    if act.ndim == 2:
        return act.float()
    if act.ndim != 3:
        return None
    B, d1, d2 = act.shape
    if name.startswith("transformer_") or name in ("trunk_proj",):
        # [B, T, D] — feature axis is last (d2)
        return act.reshape(-1, d2).float()
    else:
        # conv: [B, C, L] — feature axis is d1
        return act.permute(0, 2, 1).reshape(-1, d1).float()


def _svd_metrics(mat: torch.Tensor, top_k: int, centered: bool) -> dict:
    if centered:
        mat = mat - mat.mean(dim=0, keepdim=True)
    try:
        _, S, _ = torch.linalg.svd(mat, full_matrices=False)
    except Exception:
        return None
    S = S.numpy()
    S = S[S > 0]
    p = S / S.sum()
    effective_rank = float(np.exp(-np.sum(p * np.log(p + 1e-12))))
    energy = S ** 2
    rank_90 = int(np.searchsorted(np.cumsum(energy) / energy.sum(), 0.90) + 1)
    return {"effective_rank": effective_rank, "rank_90pct": rank_90, "singular_values": S[:top_k]}


def compute_svd_metrics(activations: Dict[str, List], top_k: int) -> List[dict]:
    results = []
    for name, tensor_list in activations.items():
        act = torch.cat(tensor_list, dim=0)
        mat = _reshape_activation(name, act)
        if mat is None:
            continue

        raw = _svd_metrics(mat, top_k, centered=False)
        cen = _svd_metrics(mat, top_k, centered=True)
        if raw is None or cen is None:
            print(f"  SVD failed for {name}")
            continue

        results.append({
            "name": name,
            "dim": mat.shape[1],
            "n_rows": mat.shape[0],
            # uncentered
            "eff_rank":    raw["effective_rank"],
            "rank_90":     raw["rank_90pct"],
            "sv":          raw["singular_values"],
            # centered (mean-subtracted per channel)
            "eff_rank_c":  cen["effective_rank"],
            "rank_90_c":   cen["rank_90pct"],
            "sv_c":        cen["singular_values"],
        })

    return results


def plot_results(results: List[dict], out_dir: str, top_k: int):
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    names       = [r["name"] for r in results]
    dims        = [r["dim"]  for r in results]
    eff_ranks   = [r["eff_rank"]   for r in results]
    ranks_90    = [r["rank_90"]    for r in results]
    eff_ranks_c = [r["eff_rank_c"] for r in results]
    ranks_90_c  = [r["rank_90_c"]  for r in results]
    x = list(range(len(names)))

    # ── Figure 1: effective rank — raw vs centered ────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(max(14, len(names) * 0.5), 9), sharex=True)
    for ax, er, r90, title in [
        (axes[0], eff_ranks,   ranks_90,   "Uncentered"),
        (axes[1], eff_ranks_c, ranks_90_c, "Centered (mean-subtracted per channel)"),
    ]:
        ax.plot(x, er,   marker="o", linewidth=1.5, label="Effective rank")
        ax.plot(x, r90,  marker="s", linewidth=1.5, linestyle="--", label="90% energy rank")
        ax.plot(x, dims, linewidth=1, linestyle=":", color="gray", label="Layer dim")
        ax.set_ylabel("Rank")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(names, rotation=90, fontsize=7)
    fig.suptitle("Per-layer effective rank (activation SVD)", y=1.01)
    fig.tight_layout()
    fig.savefig(out_path / "effective_rank.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path / 'effective_rank.png'}")

    # ── Figure 2: SV spectrum heatmap — raw vs centered ───────────────────────
    n_layers = len(results)
    k = min(top_k, min(len(r["sv"]) for r in results), min(len(r["sv_c"]) for r in results))

    fig, axes = plt.subplots(1, 2, figsize=(max(16, k * 0.2), max(6, n_layers * 0.3)))
    for ax, sv_key, title in [
        (axes[0], "sv",   "Uncentered"),
        (axes[1], "sv_c", "Centered"),
    ]:
        matrix = np.zeros((n_layers, k))
        for i, r in enumerate(results):
            sv = r[sv_key][:k]
            sv_norm = sv / (sv[0] + 1e-12)
            matrix[i, :len(sv_norm)] = sv_norm
        im = ax.imshow(matrix, aspect="auto", cmap="viridis", origin="upper")
        ax.set_yticks(range(n_layers))
        ax.set_yticklabels(names, fontsize=7)
        ax.set_xlabel("Singular value index")
        ax.set_title(f"SV spectra — {title} (σᵢ/σ₁)")
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path / "sv_spectrum.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path / 'sv_spectrum.png'}")


if __name__ == "__main__":
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"Loading model from {args.checkpoint} ...")
    model, species_configs, cfg = load_model(args.checkpoint, args.config, device)
    print(f"Model loaded on {device}")

    print(f"Loading {args.n_samples} samples for species '{args.species}' ...")
    one_hot = load_data(species_configs, args.species, cfg, args.n_samples, device)
    print(f"Data loaded: {one_hot.shape}")

    activations, handles = register_hooks(model)

    print("Running forward pass ...")
    with torch.no_grad():
        for i in range(one_hot.shape[0]):
            model(one_hot[i:i+1], head=args.species)

    for h in handles:
        h.remove()

    print(f"Collected activations for {len(activations)} layers")

    print("Computing SVD metrics ...")
    results = compute_svd_metrics(activations, args.top_k)

    header = (f"{'Layer':<35} {'Dim':>6}"
              f"  {'EffRank':>9} {'R90%':>6}"
              f"  {'EffRank(c)':>10} {'R90%(c)':>8}")
    sep = "-" * len(header)
    print("\n" + header)
    print(sep)
    for r in results:
        print(f"{r['name']:<35} {r['dim']:>6}"
              f"  {r['eff_rank']:>9.1f} {r['rank_90']:>6d}"
              f"  {r['eff_rank_c']:>10.1f} {r['rank_90_c']:>8d}")

    plot_results(results, args.out_dir, args.top_k)
    print("\nDone.")
