#!/usr/bin/env python3
"""Plot task-010 heatmaps, primary contrasts, lag structure and trajectories."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf


ALLELES = ("cpp8_1", "cpp8_3")
RT_CLASSES = ("E", "M", "L")
STAGES = ("ES", "MS", "LS")
COLORS = {"cpp8_1": "#8F3B46", "cpp8_3": "#D67B5B"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--coefficients", required=True)
    p.add_argument("--contrasts", required=True)
    p.add_argument("--lag", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir, plotdir = Path(args.output_dir), Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(args.coefficients, sep="\t")
    contrasts = pd.read_csv(args.contrasts, sep="\t")
    lag = pd.read_csv(args.lag, sep="\t")
    data = pd.read_csv(args.matrix, sep="\t")

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    vmax = float(np.nanmax(np.abs(cells["standardized_beta_delta_rt"])))
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.6), constrained_layout=True)
    for ax, allele in zip(axes, ALLELES):
        q = cells[cells["allele"] == allele]
        mat = q.pivot(index="wt_rt_class", columns="atac_stage", values="standardized_beta_delta_rt").reindex(index=RT_CLASSES, columns=STAGES)
        sns.heatmap(mat, ax=ax, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".3f", linewidths=0.5, cbar=ax is axes[-1], cbar_kws={"label": "Adjusted standardized β"})
        for ri, rt_class in enumerate(RT_CLASSES):
            for si, stage in enumerate(STAGES):
                row = q[(q["wt_rt_class"] == rt_class) & (q["atac_stage"] == stage)].iloc[0]
                if row["bootstrap_ci_low"] > 0 or row["bootstrap_ci_high"] < 0:
                    ax.add_patch(plt.Rectangle((si + 0.04, ri + 0.04), 0.92, 0.92, fill=False, edgecolor="black", lw=1.1))
        ax.set_title(allele.replace("cpp8_", "oscpp8-") + " allele")
        ax.set_xlabel("ATAC measurement stage")
        ax.set_ylabel("WT local RT class")
    fig.savefig(plotdir / "rt_class_by_atac_stage_heatmap_010.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_class_by_atac_stage_heatmap_010.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    order = ["diagonal_dominance", "post_replication_lag", "global_ms_dominance"]
    labels = ["Diagonal dominance", "+1-stage lag", "Global MS dominance"]
    for i, allele in enumerate(ALLELES):
        q = contrasts[contrasts["allele"] == allele].set_index("contrast").loc[order].reset_index()
        y = np.arange(3) + (i - 0.5) * 0.16
        axes[0].errorbar(q["estimate"], y, xerr=[q["estimate"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["estimate"]], fmt="o", capsize=3, color=COLORS[allele], label=allele.replace("cpp8_", "oscpp8-"))
    axes[0].axvline(0, color="black", lw=0.7, ls="--")
    axes[0].set_yticks(range(3), labels)
    axes[0].set_xlabel("Predefined contrast in standardized β")
    axes[0].set_title("a  Competing temporal structures", loc="left", fontweight="bold")
    axes[0].legend(frameon=False)
    for allele in ALLELES:
        q = lag[lag["allele"] == allele].sort_values("lag")
        axes[1].errorbar(q["lag"], q["mean_beta"], yerr=[q["mean_beta"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["mean_beta"]], marker="o", capsize=3, color=COLORS[allele], label=allele.replace("cpp8_", "oscpp8-"))
    axes[1].axhline(0, color="black", lw=0.7, ls="--")
    axes[1].axvline(0, color="black", lw=0.5, ls=":")
    axes[1].set_xticks(range(-2, 3))
    axes[1].set_xlabel("ATAC stage index − WT RT class index")
    axes[1].set_ylabel("Mean adjusted standardized β")
    axes[1].set_title("b  Coupling by stage distance", loc="left", fontweight="bold")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "stage_lag_contrasts_010.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "stage_lag_contrasts_010.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    trajectory_rows = []
    use = data[data["wt_rt_class"].isin(["E", "L"])].copy()
    use["shift_direction"] = np.where(use["delta_rt"] < 0, "earlier", "later")
    for keys, d in use.groupby(["allele", "wt_rt_class", "shift_direction", "atac_stage"], sort=False):
        fit = smf.ols("delta_atac ~ 1", data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
        ci = fit.conf_int().loc["Intercept"]
        trajectory_rows.append({"allele": keys[0], "wt_rt_class": keys[1], "shift_direction": keys[2], "atac_stage": keys[3], "n_ocr": len(d), "mean_delta_atac": float(fit.params["Intercept"]), "cluster_ci_low": float(ci.iloc[0]), "cluster_ci_high": float(ci.iloc[1])})
    trajectories = pd.DataFrame(trajectory_rows)
    trajectories.to_csv(outdir / "delta_atac_trajectories_010.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    for row_i, allele in enumerate(ALLELES):
        for col_i, rt_class in enumerate(("E", "L")):
            ax = axes[row_i, col_i]
            q = trajectories[(trajectories["allele"] == allele) & (trajectories["wt_rt_class"] == rt_class)]
            for direction, color in [("earlier", "#4477AA"), ("later", "#CC6677")]:
                z = q[q["shift_direction"] == direction].set_index("atac_stage").loc[list(STAGES)].reset_index()
                x = np.arange(3)
                ax.errorbar(x, z["mean_delta_atac"], yerr=[z["mean_delta_atac"] - z["cluster_ci_low"], z["cluster_ci_high"] - z["mean_delta_atac"]], marker="o", capsize=3, color=color, label=f"shifted {direction}")
            ax.axhline(0, color="black", lw=0.6, ls="--")
            ax.set_xticks(range(3), STAGES)
            ax.set_ylabel("Mean ATAC log2 fold change")
            ax.set_title(f"{allele.replace('cpp8_', 'oscpp8-')}, WT {rt_class} regions")
            ax.spines[["top", "right"]].set_visible(False)
            if row_i == 0 and col_i == 0:
                ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plotdir / "delta_atac_trajectories_010.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "delta_atac_trajectories_010.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
