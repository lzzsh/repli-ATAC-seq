#!/usr/bin/env python3
"""Plot multivariable, matched-hotspot and signed-class feature results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--models", required=True)
    p.add_argument("--matched", required=True)
    p.add_argument("--signed", required=True)
    p.add_argument("--plot-dir", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    plotdir = Path(args.plot_dir)
    plotdir.mkdir(parents=True, exist_ok=True)
    models = pd.read_csv(args.models, sep="\t")
    matched = pd.read_csv(args.matched, sep="\t")
    signed = pd.read_csv(args.signed, sep="\t")
    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, family, title in zip(axes, ["primary", "secondary_histone"], ["Primary genomic features", "Histone-mark features"]):
        q = models[models["family"] == family].sort_values("standardized_beta")
        y = np.arange(len(q))
        ax.errorbar(q["standardized_beta"], y, xerr=[q["standardized_beta"] - q["cluster_ci_low"], q["cluster_ci_high"] - q["standardized_beta"]], fmt="o", capsize=3, color="#332288")
        ax.axvline(0, color="black", lw=0.7, ls="--")
        ax.set_yticks(y, q["feature_label"])
        ax.set_xlabel("Standardized multivariable β")
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_instability_feature_forest_012.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_instability_feature_forest_012.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    q = matched[matched["hotspot_threshold_percent"] == 5].sort_values("standardized_hotspot_minus_control")
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    y = np.arange(len(q))
    ax.errorbar(q["standardized_hotspot_minus_control"], y, xerr=[q["standardized_hotspot_minus_control"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["standardized_hotspot_minus_control"]], fmt="o", capsize=3, color="#CC6677")
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.set_yticks(y, q["feature_label"])
    ax.set_xlabel("Matched hotspot − control standardized difference")
    ax.set_title("CPP8 consensus RT hotspots (top 5%)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "hotspot_matched_enrichment_012.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "hotspot_matched_enrichment_012.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    primary_labels = models.loc[models["family"] == "primary", "feature_label"].tolist()
    q = signed[(signed["feature_label"].isin(primary_labels)) & signed["shift_class_vs_stable"].isin(["concordant_earlier", "concordant_later"])].copy()
    mat = q.pivot(index="feature_label", columns="shift_class_vs_stable", values="standardized_difference").reindex(primary_labels)
    fig, ax = plt.subplots(figsize=(5.7, 5.2))
    sns.heatmap(mat, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.4, cbar_kws={"label": "Adjusted standardized difference"}, ax=ax)
    ax.set_xlabel("Consensus hotspot direction vs stable bins")
    ax.set_ylabel("")
    ax.set_title("Signed RT-shift feature profiles")
    fig.tight_layout()
    fig.savefig(plotdir / "signed_shift_feature_profiles_012.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "signed_shift_feature_profiles_012.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
