#!/usr/bin/env python3
"""Validate two-allele CPP8 paths with chromosome-LOCO prediction and spatial nulls."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


SEED = 20260826
MODELS = {
    "baseline": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"],
    "baseline_plus_rt": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_rt"],
    "baseline_plus_atac": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_atac"],
    "baseline_plus_rt_atac": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_rt", "delta_atac"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--permutations", type=int, default=500)
    return p.parse_args()


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / sd if np.isfinite(sd) and sd > 0 else np.full_like(x, np.nan)


def empirical_p(null: np.ndarray, observed: float) -> float:
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))


def qr_controls(d: pd.DataFrame, include_baseline_rna: bool) -> np.ndarray:
    cols = ["z_wt_wrt", "z_wt_log_atac", "z_log_n_ocr"]
    if include_baseline_rna:
        cols.append("z_baseline_rna")
    chrom = pd.get_dummies(d["chr"].astype(str), drop_first=True, dtype=float)
    c = np.column_stack([np.ones(len(d)), d[cols].to_numpy(float), chrom.to_numpy(float)])
    q, _ = np.linalg.qr(c, mode="reduced")
    return q


def residualize(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    return v - q @ (q.T @ v)


def simple_slope(x: np.ndarray, y: np.ndarray) -> float:
    den = float(x @ x)
    return float(x @ y / den) if den > 0 else np.nan


def mediator_slope(m: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    mm, xx, mx = float(m @ m), float(x @ x), float(m @ x)
    my, xy = float(m @ y), float(x @ y)
    den = mm * xx - mx * mx
    return float((my * xx - xy * mx) / den) if den > 0 else np.nan


def make_strata(d: pd.DataFrame, baseline: str) -> np.ndarray:
    labels = np.empty(len(d), dtype=object)
    for chrom, idx in d.groupby("chr", sort=False).groups.items():
        q = pd.qcut(d.loc[idx, baseline].rank(method="first"), 5, labels=False, duplicates="drop")
        labels[np.asarray(idx, dtype=int)] = [f"{chrom}:{int(v)}" for v in q]
    return labels.astype(str)


def stratified_permutation_matrix(values: np.ndarray, strata: np.ndarray, nperm: int, rng: np.random.Generator) -> np.ndarray:
    mat = np.empty((len(values), nperm), dtype=float)
    levels = [np.flatnonzero(strata == x) for x in np.unique(strata)]
    for i in range(nperm):
        out = values.copy()
        for idx in levels:
            if len(idx) > 1:
                out[idx] = values[rng.permutation(idx)]
        mat[:, i] = out
    return mat


def circular_permutation_matrix(values: np.ndarray, d: pd.DataFrame, nperm: int, rng: np.random.Generator) -> np.ndarray:
    mat = np.empty((len(values), nperm), dtype=float)
    ordered_by_chr = []
    for _, idx in d.groupby("chr", sort=False).groups.items():
        ordered_by_chr.append(np.asarray(sorted(idx, key=lambda j: d.at[j, "tss"]), dtype=int))
    for i in range(nperm):
        out = values.copy()
        for ordered in ordered_by_chr:
            n = len(ordered)
            if n < 3:
                continue
            low = max(1, n // 10)
            high = max(low + 1, n - low)
            out[ordered] = np.roll(values[ordered], int(rng.integers(low, high)))
        mat[:, i] = out
    return mat


def path_nulls(d: pd.DataFrame, nperm: int, rng: np.random.Generator) -> list[dict[str, object]]:
    d = d.reset_index(drop=True).copy()
    cols = ["delta_rt", "delta_atac", "rna_logfc", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]
    for col in cols:
        d[f"z_{col}"] = zscore(d[col].to_numpy(float))
    x = d["z_delta_rt"].to_numpy(float)
    m = d["z_delta_atac"].to_numpy(float)
    y = d["z_rna_logfc"].to_numpy(float)
    qa = qr_controls(d, include_baseline_rna=False)
    qb = qr_controls(d, include_baseline_rna=True)
    rx_a, rm_a = residualize(x, qa), residualize(m, qa)
    rx_b, rm_b, ry_b = residualize(x, qb), residualize(m, qb), residualize(y, qb)
    observed = simple_slope(rx_a, rm_a) * mediator_slope(rm_b, rx_b, ry_b)

    def eval_rt(xp: np.ndarray) -> np.ndarray:
        rxpa = xp - qa @ (qa.T @ xp)
        rxpb = xp - qb @ (qb.T @ xp)
        a = np.sum(rxpa * rm_a[:, None], axis=0) / np.sum(rxpa * rxpa, axis=0)
        mm, my = float(rm_b @ rm_b), float(rm_b @ ry_b)
        xx = np.sum(rxpb * rxpb, axis=0)
        mx = np.sum(rm_b[:, None] * rxpb, axis=0)
        xy = np.sum(ry_b[:, None] * rxpb, axis=0)
        b = (my * xx - xy * mx) / (mm * xx - mx * mx)
        return a * b

    def eval_atac(mp: np.ndarray) -> np.ndarray:
        rma = mp - qa @ (qa.T @ mp)
        rmb = mp - qb @ (qb.T @ mp)
        a = np.sum(rx_a[:, None] * rma, axis=0) / float(rx_a @ rx_a)
        mm = np.sum(rmb * rmb, axis=0)
        xx, xy = float(rx_b @ rx_b), float(rx_b @ ry_b)
        mx = np.sum(rmb * rx_b[:, None], axis=0)
        my = np.sum(rmb * ry_b[:, None], axis=0)
        b = (my * xx - xy * mx) / (mm * xx - mx * mx)
        return a * b

    rt_strata = make_strata(d, "wt_wrt")
    atac_strata = make_strata(d, "wt_log_atac")
    nulls = {
        "rt_within_chr_wrt_quintile": eval_rt(stratified_permutation_matrix(x, rt_strata, nperm, rng)),
        "rt_chr_circular_shift": eval_rt(circular_permutation_matrix(x, d, nperm, rng)),
        "atac_within_chr_atac_quintile": eval_atac(stratified_permutation_matrix(m, atac_strata, nperm, rng)),
        "atac_chr_circular_shift": eval_atac(circular_permutation_matrix(m, d, nperm, rng)),
    }
    rows = []
    for scheme, null in nulls.items():
        rows.append({
            "scheme": scheme, "n_permutations": nperm, "observed_indirect": observed,
            "null_mean": float(np.mean(null)), "null_sd": float(np.std(null, ddof=1)),
            "null_ci_low": float(np.percentile(null, 2.5)), "null_ci_high": float(np.percentile(null, 97.5)),
            "empirical_p_two_sided": empirical_p(null, observed),
        })
    return rows


def chromosome_loco(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    folds = []
    predictions = []
    for model_name, predictors in MODELS.items():
        pred = np.full(len(d), np.nan)
        for chrom in sorted(d["chr"].unique()):
            test = d["chr"] == chrom
            train = ~test
            model = make_pipeline(StandardScaler(), LinearRegression())
            model.fit(d.loc[train, predictors], d.loc[train, "rna_logfc"])
            pred[test] = model.predict(d.loc[test, predictors])
            folds.append({"model": model_name, "left_out_chromosome": chrom, "n_test": int(test.sum()), "r2": r2_score(d.loc[test, "rna_logfc"], pred[test]), "mse": mean_squared_error(d.loc[test, "rna_logfc"], pred[test])})
        predictions.append(pd.DataFrame({"model": model_name, "observed": d["rna_logfc"], "predicted": pred}))
    folds = pd.DataFrame(folds)
    pred = pd.concat(predictions, ignore_index=True)
    summary = []
    for model_name in MODELS:
        q = pred[pred["model"] == model_name]
        summary.append({"model": model_name, "pooled_loco_r2": r2_score(q["observed"], q["predicted"]), "pooled_loco_rmse": float(np.sqrt(mean_squared_error(q["observed"], q["predicted"]))), "median_chromosome_r2": float(folds.loc[folds["model"] == model_name, "r2"].median())})
    summary = pd.DataFrame(summary)
    baseline = float(summary.loc[summary["model"] == "baseline", "pooled_loco_r2"].iloc[0])
    summary["delta_pooled_r2_vs_baseline"] = summary["pooled_loco_r2"] - baseline
    comparisons = []
    for a, b, label in [
        ("baseline", "baseline_plus_rt", "add_rt_to_baseline"),
        ("baseline", "baseline_plus_atac", "add_atac_to_baseline"),
        ("baseline", "baseline_plus_rt_atac", "add_rt_and_atac_to_baseline"),
        ("baseline_plus_atac", "baseline_plus_rt_atac", "add_rt_given_atac"),
        ("baseline_plus_rt", "baseline_plus_rt_atac", "add_atac_given_rt"),
    ]:
        ma = folds[folds["model"] == a].set_index("left_out_chromosome")["mse"]
        mb = folds[folds["model"] == b].set_index("left_out_chromosome")["mse"]
        diff = ma - mb
        try:
            p = float(wilcoxon(diff, alternative="greater").pvalue)
        except ValueError:
            p = np.nan
        comparisons.append({"comparison": label, "mean_mse_reduction": float(diff.mean()), "chromosome_wilcoxon_p": p})
    return summary, folds, pd.DataFrame(comparisons)


def main() -> None:
    args = parse_args()
    outdir, plotdir = Path(args.output_dir), Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.gene_metrics, sep="\t")
    metrics["log_n_ocr"] = np.log1p(metrics["n_promoter_ocr"])

    summaries, folds_all, comparisons_all, null_rows = [], [], [], []
    for keys, d in metrics.groupby(["allele", "stage", "promoter_window_bp"], sort=False):
        allele, stage, window = keys
        d = d.dropna(subset=["rna_logfc", "delta_rt", "delta_atac", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]).reset_index(drop=True)
        summary, folds, comparisons = chromosome_loco(d)
        for frame in [summary, folds, comparisons]:
            frame["allele"], frame["stage"], frame["promoter_window_bp"] = allele, stage, int(window)
        summaries.append(summary)
        folds_all.append(folds)
        comparisons_all.append(comparisons)
        for row in path_nulls(d, args.permutations, rng):
            row.update({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "n_genes": len(d)})
            null_rows.append(row)

    summary = pd.concat(summaries, ignore_index=True)
    folds = pd.concat(folds_all, ignore_index=True)
    comparisons = pd.concat(comparisons_all, ignore_index=True)
    nulls = pd.DataFrame(null_rows)
    nulls["fdr_within_permutation_scheme"] = nulls.groupby("scheme")["empirical_p_two_sided"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    summary.to_csv(outdir / "cpp8_allele_loco_009.tsv", sep="\t", index=False)
    folds.to_csv(outdir / "cpp8_allele_loco_folds_009.tsv", sep="\t", index=False)
    comparisons.to_csv(outdir / "cpp8_allele_loco_comparisons_009.tsv", sep="\t", index=False)
    nulls.to_csv(outdir / "cpp8_allele_spatial_null_009.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    plot_cv = summary[summary["model"] != "baseline"].copy()
    plot_cv["label"] = plot_cv["allele"].str.replace("cpp8_", "oscpp8-") + " " + plot_cv["stage"] + " ±" + (plot_cv["promoter_window_bp"] // 1000).astype(str) + "kb"
    labels = list(dict.fromkeys(plot_cv["label"]))
    ypos = {label: i for i, label in enumerate(labels)}
    offsets = {"baseline_plus_rt": -0.18, "baseline_plus_atac": 0.0, "baseline_plus_rt_atac": 0.18}
    colors = {"baseline_plus_rt": "#AA4499", "baseline_plus_atac": "#44AA99", "baseline_plus_rt_atac": "#4477AA"}
    for model, q in plot_cv.groupby("model"):
        axes[0].scatter(q["delta_pooled_r2_vs_baseline"], [ypos[x] + offsets[model] for x in q["label"]], s=18, color=colors[model], label=model.replace("baseline_plus_", "+"))
    axes[0].axvline(0, color="black", lw=0.7, ls="--")
    axes[0].set_yticks(range(len(labels)), labels)
    axes[0].set_xlabel("Change in chromosome-LOCO R² vs baseline")
    axes[0].set_title("a  Out-of-chromosome predictive contribution", loc="left", fontweight="bold")
    axes[0].legend(frameon=False)

    q = nulls[nulls["scheme"].isin(["rt_chr_circular_shift", "atac_chr_circular_shift"])].copy()
    q["label"] = q["allele"].str.replace("cpp8_", "oscpp8-") + " " + q["stage"] + " ±" + (q["promoter_window_bp"] // 1000).astype(str) + "kb"
    for i, scheme in enumerate(["rt_chr_circular_shift", "atac_chr_circular_shift"]):
        z = q[q["scheme"] == scheme]
        axes[1].scatter(-np.log10(z["empirical_p_two_sided"]), [ypos[x] + (-0.08 if i == 0 else 0.08) for x in z["label"]], s=18, label=scheme.replace("_chr_circular_shift", " shift"))
    axes[1].axvline(-np.log10(0.05), color="black", lw=0.7, ls="--")
    axes[1].set_yticks(range(len(labels)), labels)
    axes[1].set_xlabel("−log10 empirical P")
    axes[1].set_title("b  Spatial circular-shift null tests", loc="left", fontweight="bold")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "cpp8_allele_validation_009.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "cpp8_allele_validation_009.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("LOCO summaries")
    print(summary.to_string(index=False))
    print("\nSpatial nulls")
    print(nulls.to_string(index=False))


if __name__ == "__main__":
    main()
