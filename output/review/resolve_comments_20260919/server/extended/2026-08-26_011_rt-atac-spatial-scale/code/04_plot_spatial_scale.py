#!/usr/bin/env python3
"""Plot nested-scale decay, offset curves and RT-boundary sensitivity."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ALLELES = ("cpp8_1", "cpp8_3")
STAGES = ("ES", "MS", "LS")
COLORS = {"ES": "#4477AA", "MS": "#228833", "LS": "#CC6677"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--coefficients", required=True)
    p.add_argument("--offsets", required=True)
    p.add_argument("--boundaries", required=True)
    p.add_argument("--plot-dir", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    plotdir = Path(args.plot_dir)
    plotdir.mkdir(parents=True, exist_ok=True)
    coefficients = pd.read_csv(args.coefficients, sep="\t")
    offsets = pd.read_csv(args.offsets, sep="\t")
    boundaries = pd.read_csv(args.boundaries, sep="\t")
    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    for ax, allele in zip(axes, ALLELES):
        for stage in STAGES:
            q = coefficients[(coefficients["allele"] == allele) & (coefficients["stage"] == stage)].sort_values("scale_radius_kb")
            ax.errorbar(q["scale_radius_kb"], q["standardized_beta_delta_rt"], yerr=[q["standardized_beta_delta_rt"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["standardized_beta_delta_rt"]], marker="o", capsize=3, color=COLORS[stage], label=stage)
        ax.axhline(0, color="black", lw=0.7, ls="--")
        ax.set_xscale("log")
        ax.set_xticks([1, 5, 25, 100], ["1", "5", "25", "100"])
        ax.set_xlabel("RT averaging radius (kb)")
        ax.set_title(allele.replace("cpp8_", "oscpp8-") + " allele")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Adjusted standardized β: ΔRT → ΔATAC")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_atac_scale_decay_011.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_atac_scale_decay_011.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, allele in zip(axes, ALLELES):
        for stage in STAGES:
            q = offsets[(offsets["allele"] == allele) & (offsets["stage"] == stage)].sort_values("offset_kb")
            ax.plot(q["offset_kb"], q["standardized_beta_delta_rt"], color=COLORS[stage], label=stage)
            ax.fill_between(q["offset_kb"], q["cluster_ci_low"], q["cluster_ci_high"], color=COLORS[stage], alpha=0.12)
        ax.axhline(0, color="black", lw=0.7, ls="--")
        ax.axvline(0, color="black", lw=0.7, ls=":")
        ax.set_xlabel("RT window offset from OCR midpoint (kb)")
        ax.set_title(allele.replace("cpp8_", "oscpp8-") + " allele")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Adjusted standardized β")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_atac_offset_curve_011.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_atac_offset_curve_011.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    primary = boundaries[(boundaries["boundary_top_percent"] == 10) & boundaries["effect"].isin(["interior", "near", "near_minus_interior"])].copy()
    labels = [f"{a.replace('cpp8_', 'oscpp8-')} {s}" for a in ALLELES for s in STAGES]
    ypos = {label: i for i, label in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for effect, marker, color, offset in [("interior", "o", "#999999", -0.18), ("near", "o", "#332288", 0.0), ("near_minus_interior", "D", "#CC6677", 0.18)]:
        q = primary[primary["effect"] == effect]
        q = q.assign(label=q["allele"].str.replace("cpp8_", "oscpp8-") + " " + q["stage"])
        y = [ypos[x] + offset for x in q["label"]]
        ax.errorbar(q["estimate"], y, xerr=[q["estimate"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["estimate"]], fmt=marker, capsize=3, color=color, label=effect.replace("_", " "))
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("Adjusted standardized β or near−interior difference")
    ax.set_title("Candidate RT-boundary coupling (top 10% local maxima)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_boundary_coupling_011.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_boundary_coupling_011.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
