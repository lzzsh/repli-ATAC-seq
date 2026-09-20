#!/usr/bin/env python3
"""Validate RT-ATAC-RNA path results with out-of-chromosome prediction and null tests."""

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
MODEL_COLORS = {
    "baseline_plus_rt": "#AA4499",
    "baseline_plus_atac": "#44AA99",
    "baseline_plus_rt_atac": "#4477AA",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--path-effects", required=True)
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
    mm = float(m @ m)
    xx = float(x @ x)
    mx = float(m @ x)
    my = float(m @ y)
    xy = float(x @ y)
    den = mm * xx - mx * mx
    return float((my * xx - xy * mx) / den) if den > 0 else np.nan


def stratified_permute(values: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = values.copy()
    for level in np.unique(strata):
        idx = np.flatnonzero(strata == level)
        if len(idx) > 1:
            out[idx] = values[rng.permutation(idx)]
    return out


def circular_permute(values: np.ndarray, d: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    out = values.copy()
    for _, idx in d.groupby("chr", sort=False).groups.items():
        ordered = np.asarray(sorted(idx, key=lambda i: d.at[i, "tss"]), dtype=int)
        n = len(ordered)
        if n < 3:
            continue
        low = max(1, n // 10)
        high = max(low + 1, n - low)
        shift = int(rng.integers(low, high))
        out[ordered] = np.roll(values[ordered], shift)
    return out


def make_strata(d: pd.DataFrame, baseline: str) -> np.ndarray:
    labels = np.empty(len(d), dtype=object)
    for chrom, idx in d.groupby("chr", sort=False).groups.items():
        x = d.loc[idx, baseline]
        q = pd.qcut(x.rank(method="first"), 5, labels=False, duplicates="drop")
        labels[np.asarray(idx, dtype=int)] = [f"{chrom}:{int(v)}" for v in q]
    return labels.astype(str)


def path_permutation_tests(d: pd.DataFrame, nperm: int, rng: np.random.Generator) -> list[dict[str, object]]:
    d = d.reset_index(drop=True).copy()
    for col in ["delta_rt", "delta_atac", "rna_logfc", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]:
        d[f"z_{col}"] = zscore(d[col].to_numpy(float))
    x = d["z_delta_rt"].to_numpy(float)
    m = d["z_delta_atac"].to_numpy(float)
    y = d["z_rna_logfc"].to_numpy(float)
    qa = qr_controls(d, include_baseline_rna=False)
    qb = qr_controls(d, include_baseline_rna=True)
    rx_a = residualize(x, qa)
    rm_a = residualize(m, qa)
    rx_b = residualize(x, qb)
    rm_b = residualize(m, qb)
    ry_b = residualize(y, qb)
    observed = simple_slope(rx_a, rm_a) * mediator_slope(rm_b, rx_b, ry_b)
    rt_strata = make_strata(d, "wt_wrt")
    atac_strata = make_strata(d, "wt_log_atac")
    def permutation_matrix(values: np.ndarray, mode: str, strata: np.ndarray | None = None) -> np.ndarray:
        mat = np.empty((len(values), nperm), dtype=float)
        for i in range(nperm):
            if mode == "stratified":
                assert strata is not None
                mat[:, i] = stratified_permute(values, strata, rng)
            else:
                mat[:, i] = circular_permute(values, d, rng)
        return mat

    def evaluate_rt_permutations(xp: np.ndarray) -> np.ndarray:
        rxp_a = xp - qa @ (qa.T @ xp)
        rxp_b = xp - qb @ (qb.T @ xp)
        a = np.sum(rxp_a * rm_a[:, None], axis=0) / np.sum(rxp_a * rxp_a, axis=0)
        mm = float(rm_b @ rm_b)
        my = float(rm_b @ ry_b)
        xx = np.sum(rxp_b * rxp_b, axis=0)
        mx = np.sum(rm_b[:, None] * rxp_b, axis=0)
        xy = np.sum(ry_b[:, None] * rxp_b, axis=0)
        b = (my * xx - xy * mx) / (mm * xx - mx * mx)
        return a * b

    def evaluate_atac_permutations(mp: np.ndarray) -> np.ndarray:
        rmp_a = mp - qa @ (qa.T @ mp)
        rmp_b = mp - qb @ (qb.T @ mp)
        a = np.sum(rx_a[:, None] * rmp_a, axis=0) / float(rx_a @ rx_a)
        mm = np.sum(rmp_b * rmp_b, axis=0)
        xx = float(rx_b @ rx_b)
        mx = np.sum(rmp_b * rx_b[:, None], axis=0)
        my = np.sum(rmp_b * ry_b[:, None], axis=0)
        xy = float(rx_b @ ry_b)
        b = (my * xx - xy * mx) / (mm * xx - mx * mx)
        return a * b


    nulls = {}
    nulls["rt_within_chr_wrt_quintile"] = evaluate_rt_permutations(
        permutation_matrix(x, "stratified", rt_strata)
    )
    nulls["rt_chr_circular_shift"] = evaluate_rt_permutations(
        permutation_matrix(x, "circular")
    )
    nulls["atac_within_chr_atac_quintile"] = evaluate_atac_permutations(
        permutation_matrix(m, "stratified", atac_strata)
    )
    nulls["atac_chr_circular_shift"] = evaluate_atac_permutations(
        permutation_matrix(m, "circular")
    )
    rows = []
    for scheme, null in nulls.items():
        rows.append(
            {
                "scheme": scheme,
                "n_permutations": nperm,
                "observed_indirect": observed,
                "null_mean": float(np.mean(null)),
                "null_sd": float(np.std(null, ddof=1)),
                "null_ci_low": float(np.percentile(null, 2.5)),
                "null_ci_high": float(np.percentile(null, 97.5)),
                "empirical_p_two_sided": empirical_p(null, observed),
            }
        )
    return rows


def chromosome_loco_cv(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fold_rows = []
    pred_frames = []
    for model_name, predictors in MODELS.items():
        pred = np.full(len(d), np.nan)
        for chrom in sorted(d["chr"].unique()):
            test = d["chr"] == chrom
            train = ~test
            model = make_pipeline(StandardScaler(), LinearRegression())
            model.fit(d.loc[train, predictors], d.loc[train, "rna_logfc"])
            pred[test] = model.predict(d.loc[test, predictors])
            fold_rows.append(
                {
                    "model": model_name,
                    "left_out_chromosome": chrom,
                    "n_test": int(test.sum()),
                    "r2": r2_score(d.loc[test, "rna_logfc"], pred[test]),
                    "mse": mean_squared_error(d.loc[test, "rna_logfc"], pred[test]),
                }
            )
        pred_frames.append(pd.DataFrame({"gene_id": d["gene_id"], "chr": d["chr"], "model": model_name, "observed": d["rna_logfc"], "predicted": pred}))
    folds = pd.DataFrame(fold_rows)
    predictions = pd.concat(pred_frames, ignore_index=True)
    summary_rows = []
    for model_name in MODELS:
        q = predictions[predictions["model"] == model_name]
        summary_rows.append(
            {
                "model": model_name,
                "pooled_loco_r2": r2_score(q["observed"], q["predicted"]),
                "pooled_loco_rmse": float(np.sqrt(mean_squared_error(q["observed"], q["predicted"]))),
                "median_chromosome_r2": float(folds.loc[folds.model == model_name, "r2"].median()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    base_r2 = float(summary.loc[summary.model == "baseline", "pooled_loco_r2"].iloc[0])
    summary["delta_pooled_r2_vs_baseline"] = summary["pooled_loco_r2"] - base_r2
    comparisons = []
    for a, b, label in [
        ("baseline", "baseline_plus_rt", "add_rt_to_baseline"),
        ("baseline", "baseline_plus_atac", "add_atac_to_baseline"),
        ("baseline", "baseline_plus_rt_atac", "add_rt_and_atac_to_baseline"),
        ("baseline_plus_atac", "baseline_plus_rt_atac", "add_rt_given_atac"),
        ("baseline_plus_rt", "baseline_plus_rt_atac", "add_atac_given_rt"),
    ]:
        ma = folds[folds.model == a].set_index("left_out_chromosome")["mse"]
        mb = folds[folds.model == b].set_index("left_out_chromosome")["mse"]
        diff = ma - mb
        try:
            pvalue = float(wilcoxon(diff, alternative="greater").pvalue)
        except ValueError:
            pvalue = np.nan
        comparisons.append({"comparison": label, "mean_mse_reduction": float(diff.mean()), "chromosome_wilcoxon_p": pvalue})
    return summary, folds, pd.DataFrame(comparisons)


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.gene_metrics, sep="\t")
    metrics["log_n_ocr"] = np.log1p(metrics["n_promoter_ocr"])
    effects = pd.read_csv(args.path_effects, sep="\t")

    cv_summary_frames = []
    cv_fold_frames = []
    cv_compare_frames = []
    permutation_rows = []
    group_cols = ["mutant", "stage", "promoter_window_bp"]
    for keys, d in metrics.groupby(group_cols, sort=False):
        mutant, stage, window = keys
        d = d.dropna(subset=["rna_logfc", "delta_rt", "delta_atac", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]).reset_index(drop=True)
        summary, folds, comparisons = chromosome_loco_cv(d)
        for frame in [summary, folds, comparisons]:
            frame["mutant"] = mutant
            frame["stage"] = stage
            frame["promoter_window_bp"] = int(window)
        cv_summary_frames.append(summary)
        cv_fold_frames.append(folds)
        cv_compare_frames.append(comparisons)
        for row in path_permutation_tests(d, args.permutations, rng):
            row.update({"mutant": mutant, "stage": stage, "promoter_window_bp": int(window), "n_genes": len(d)})
            permutation_rows.append(row)

    cv_summary = pd.concat(cv_summary_frames, ignore_index=True)
    cv_folds = pd.concat(cv_fold_frames, ignore_index=True)
    cv_comparisons = pd.concat(cv_compare_frames, ignore_index=True)
    perm = pd.DataFrame(permutation_rows)
    perm["fdr_within_permutation_scheme"] = perm.groupby("scheme")["empirical_p_two_sided"].transform(
        lambda x: multipletests(x, method="fdr_bh")[1]
    )
    cv_summary.to_csv(outdir / "chromosome_loco_cv_summary_006.tsv", sep="\t", index=False)
    cv_folds.to_csv(outdir / "chromosome_loco_cv_folds_006.tsv", sep="\t", index=False)
    cv_comparisons.to_csv(outdir / "chromosome_loco_cv_model_comparisons_006.tsv", sep="\t", index=False)
    perm.to_csv(outdir / "causal_path_permutation_tests_006.tsv", sep="\t", index=False)

    wide = effects.pivot_table(
        index=["mutant", "stage"], columns="promoter_window_bp",
        values=["indirect_effect_a_times_b", "indirect_boot_ci_low", "indirect_boot_ci_high", "indirect_fdr_across_12_tests"],
    )
    wide.columns = [f"{metric}_{int(window)}bp" for metric, window in wide.columns]
    wide = wide.reset_index()
    wide["same_effect_sign_1kb_3kb"] = np.sign(wide["indirect_effect_a_times_b_1000bp"]) == np.sign(wide["indirect_effect_a_times_b_3000bp"])
    wide["both_bootstrap_ci_exclude_zero"] = (
        ((wide["indirect_boot_ci_low_1000bp"] > 0) | (wide["indirect_boot_ci_high_1000bp"] < 0))
        & ((wide["indirect_boot_ci_low_3000bp"] > 0) | (wide["indirect_boot_ci_high_3000bp"] < 0))
    )
    wide.to_csv(outdir / "promoter_window_sensitivity_006.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    plot_cv = cv_summary[cv_summary.model != "baseline"].copy()
    plot_cv["label"] = plot_cv["mutant"].str.upper() + " " + plot_cv["stage"] + " ±" + (plot_cv["promoter_window_bp"] // 1000).astype(str) + "kb"
    labels = list(dict.fromkeys(plot_cv["label"]))
    ypos = {label: i for i, label in enumerate(labels)}
    offsets = {"baseline_plus_rt": -0.18, "baseline_plus_atac": 0.0, "baseline_plus_rt_atac": 0.18}
    for model, md in plot_cv.groupby("model"):
        axes[0].scatter(md["delta_pooled_r2_vs_baseline"], [ypos[x] + offsets[model] for x in md["label"]], s=18, color=MODEL_COLORS[model], label=model.replace("baseline_plus_", "+"))
    axes[0].axvline(0, color="black", lw=0.7, ls="--")
    axes[0].set_yticks(range(len(labels)), labels)
    axes[0].set_xlabel("Change in chromosome-LOCO R² vs baseline")
    axes[0].set_title("a  Out-of-chromosome predictive contribution", loc="left", fontweight="bold")
    axes[0].legend(frameon=False)

    primary = perm[perm["scheme"].isin(["rt_chr_circular_shift", "atac_chr_circular_shift"])].copy()
    primary["label"] = primary["mutant"].str.upper() + " " + primary["stage"] + " ±" + (primary["promoter_window_bp"] // 1000).astype(str) + "kb"
    schemes = ["rt_chr_circular_shift", "atac_chr_circular_shift"]
    for i, scheme in enumerate(schemes):
        q = primary[primary.scheme == scheme]
        axes[1].scatter(-np.log10(q["empirical_p_two_sided"]), [ypos[x] + (-0.08 if i == 0 else 0.08) for x in q["label"]], s=18, label=scheme.replace("_chr_circular_shift", " shift"))
    axes[1].axvline(-np.log10(0.05), color="black", lw=0.7, ls="--")
    axes[1].set_yticks(range(len(labels)), labels)
    axes[1].set_xlabel("−log10 empirical P")
    axes[1].set_title("b  Spatial circular-shift null tests", loc="left", fontweight="bold")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "causal_path_validation_006.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "causal_path_validation_006.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Chromosome-LOCO CV")
    print(cv_summary.to_string(index=False))
    print("\nPermutation tests")
    print(perm.to_string(index=False))


if __name__ == "__main__":
    main()
