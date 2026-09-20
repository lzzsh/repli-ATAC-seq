#!/usr/bin/env python3
"""Quantify cross-dataset reproducibility of Repli-ATAC RT assignments."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score


DATASETS = ["NIP.1", "ZH11.1", "ZH11.2"]
PHASE_COLS = [f"{x}_phase" for x in DATASETS]
STAGES = ("ES", "MS", "LS")
FINAL_TO_CALL = {
    "E": "ES",
    "EM": "ESMS",
    "M": "MS",
    "ML": "MSLS",
    "L": "LS",
    "EL": "ESLS",
    "EML": "ESMSLS",
}
FINAL_ORDER = ["E", "EM", "M", "ML", "L", "EL", "EML"]
COLORS = {
    "E": "#2C5F9E",
    "EM": "#68A0D8",
    "M": "#95BE6C",
    "ML": "#E4B660",
    "L": "#E68364",
    "EL": "#9A4D8E",
    "EML": "#767676",
}


def phase_set(label: object) -> set[str]:
    if not isinstance(label, str) or label in {"Non-replication", "unknown", "nan"}:
        return set()
    return {stage for stage in STAGES if stage in label}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plot-dir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)

    signal_cols = [f"{dataset}.{stage}_norm" for dataset in DATASETS for stage in STAGES]
    usecols = ["chr", "start", "end", "final_phase", *PHASE_COLS, *signal_cols]
    df = pd.read_csv(args.input, usecols=usecols)
    df["start"] = pd.to_numeric(df["start"].astype(str).str.strip(), errors="coerce")
    df["end"] = pd.to_numeric(df["end"].astype(str).str.strip(), errors="coerce")

    pair_rows: list[dict[str, object]] = []
    for a, b in combinations(DATASETS, 2):
        ca, cb = f"{a}_phase", f"{b}_phase"
        valid = df[ca].notna() & df[cb].notna()
        labels_a = df.loc[valid, ca].astype(str)
        labels_b = df.loc[valid, cb].astype(str)
        js = [jaccard(phase_set(x), phase_set(y)) for x, y in zip(labels_a, labels_b)]
        pair_rows.append(
            {
                "dataset_a": a,
                "dataset_b": b,
                "n": int(valid.sum()),
                "exact_agreement": float((labels_a.values == labels_b.values).mean()),
                "cohen_kappa": float(cohen_kappa_score(labels_a, labels_b)),
                "mean_stage_jaccard": float(np.mean(js)),
            }
        )
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(outdir / "rt_map_pairwise_agreement_002.tsv", sep="\t", index=False)

    signal_rows: list[dict[str, object]] = []
    for stage in STAGES:
        for a, b in combinations(DATASETS, 2):
            xa = pd.to_numeric(df[f"{a}.{stage}_norm"], errors="coerce")
            xb = pd.to_numeric(df[f"{b}.{stage}_norm"], errors="coerce")
            valid = xa.notna() & xb.notna()
            rho, pvalue = spearmanr(xa[valid], xb[valid])
            signal_rows.append(
                {
                    "stage": stage,
                    "dataset_a": a,
                    "dataset_b": b,
                    "n": int(valid.sum()),
                    "spearman_rho": float(rho),
                    "p_value": float(pvalue),
                }
            )
    signal_df = pd.DataFrame(signal_rows)
    signal_df.to_csv(outdir / "rt_map_signal_correlations_002.tsv", sep="\t", index=False)

    calls = df[PHASE_COLS].astype(str).to_numpy()
    final_calls = df["final_phase"].map(FINAL_TO_CALL).fillna("").to_numpy()
    unanimous = np.array([len(set(row)) == 1 for row in calls])
    max_identical = np.array([max(pd.Series(row).value_counts().values) for row in calls])
    exact_support = np.array([(row == final).sum() for row, final in zip(calls, final_calls)])

    mean_jaccard = []
    stage_vote_min = []
    for row, final in zip(calls, final_calls):
        target = phase_set(final)
        call_sets = [phase_set(x) for x in row]
        mean_jaccard.append(np.mean([jaccard(target, x) for x in call_sets]))
        if target:
            votes = [sum(stage in x for x in call_sets) for stage in target]
            stage_vote_min.append(min(votes))
        else:
            stage_vote_min.append(0)

    confidence = df[["chr", "start", "end", "final_phase", *PHASE_COLS]].copy()
    confidence["unanimous_exact"] = unanimous
    confidence["max_identical_calls"] = max_identical
    confidence["exact_final_support"] = exact_support
    confidence["minimum_final_stage_votes"] = stage_vote_min
    confidence["mean_stage_jaccard_to_final"] = mean_jaccard
    confidence.to_csv(outdir / "rt_map_region_confidence_002.tsv.gz", sep="\t", index=False)

    by_class = (
        confidence.groupby("final_phase", observed=False)
        .agg(
            n=("final_phase", "size"),
            unanimous_rate=("unanimous_exact", "mean"),
            exact_majority_rate=("max_identical_calls", lambda x: np.mean(np.asarray(x) >= 2)),
            mean_exact_final_support=("exact_final_support", "mean"),
            mean_minimum_stage_votes=("minimum_final_stage_votes", "mean"),
            mean_stage_jaccard=("mean_stage_jaccard_to_final", "mean"),
        )
        .reset_index()
    )
    by_class["final_phase"] = pd.Categorical(by_class["final_phase"], FINAL_ORDER, ordered=True)
    by_class = by_class.sort_values("final_phase")
    by_class.to_csv(outdir / "rt_map_confidence_by_class_002.tsv", sep="\t", index=False)

    overall = pd.DataFrame(
        [
            {"metric": "regions", "value": len(df)},
            {"metric": "unanimous_exact_fraction", "value": unanimous.mean()},
            {"metric": "two_or_more_identical_fraction", "value": np.mean(max_identical >= 2)},
            {"metric": "no_exact_majority_fraction", "value": np.mean(max_identical < 2)},
            {"metric": "mean_stage_jaccard_to_final", "value": np.mean(mean_jaccard)},
            {"metric": "final_components_supported_by_all_three_fraction", "value": np.mean(np.asarray(stage_vote_min) == 3)},
        ]
    )
    overall.to_csv(outdir / "rt_map_reproducibility_summary_002.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5))

    counts = df["final_phase"].value_counts().reindex(FINAL_ORDER).dropna()
    axes[0, 0].bar(counts.index, counts.values, color=[COLORS[x] for x in counts.index])
    axes[0, 0].set_title("a  Final RT-class distribution", loc="left", fontweight="bold")
    axes[0, 0].set_ylabel("OCR count")
    axes[0, 0].tick_params(axis="x", rotation=45)

    kappa = np.eye(len(DATASETS))
    for row in pair_rows:
        i, j = DATASETS.index(row["dataset_a"]), DATASETS.index(row["dataset_b"])
        kappa[i, j] = kappa[j, i] = row["cohen_kappa"]
    im = axes[0, 1].imshow(kappa, vmin=0, vmax=1, cmap="Blues")
    axes[0, 1].set_xticks(range(len(DATASETS)), DATASETS, rotation=35, ha="right")
    axes[0, 1].set_yticks(range(len(DATASETS)), DATASETS)
    for i in range(len(DATASETS)):
        for j in range(len(DATASETS)):
            axes[0, 1].text(j, i, f"{kappa[i, j]:.3f}", ha="center", va="center")
    axes[0, 1].set_title("b  Pairwise phase-call Cohen's κ", loc="left", fontweight="bold")
    fig.colorbar(im, ax=axes[0, 1], fraction=0.046, pad=0.04)

    support_labels = ["3 identical", "2 identical", "No exact majority"]
    support_values = [np.mean(max_identical == 3), np.mean(max_identical == 2), np.mean(max_identical == 1)]
    axes[1, 0].bar(support_labels, np.asarray(support_values) * 100, color=["#2C5F9E", "#95BE6C", "#E68364"])
    axes[1, 0].set_ylabel("Regions (%)")
    axes[1, 0].set_ylim(0, 100)
    axes[1, 0].set_title("c  Agreement among three classified datasets", loc="left", fontweight="bold")
    axes[1, 0].tick_params(axis="x", rotation=20)

    plot_data = [
        confidence.loc[confidence["final_phase"] == phase, "mean_stage_jaccard_to_final"].to_numpy()
        for phase in FINAL_ORDER
        if (confidence["final_phase"] == phase).any()
    ]
    plot_labels = [phase for phase in FINAL_ORDER if (confidence["final_phase"] == phase).any()]
    bp = axes[1, 1].boxplot(plot_data, labels=plot_labels, showfliers=False, patch_artist=True)
    for patch, phase in zip(bp["boxes"], plot_labels):
        patch.set_facecolor(COLORS[phase])
        patch.set_alpha(0.8)
    axes[1, 1].set_ylim(0, 1.03)
    axes[1, 1].set_ylabel("Mean stage-set Jaccard to final call")
    axes[1, 1].set_title("d  Region-level RT confidence", loc="left", fontweight="bold")

    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_map_reproducibility_002.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_map_reproducibility_002.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(overall.to_string(index=False))
    print("\nPairwise agreement")
    print(pair_df.to_string(index=False))
    print("\nBy class")
    print(by_class.to_string(index=False))


if __name__ == "__main__":
    main()
