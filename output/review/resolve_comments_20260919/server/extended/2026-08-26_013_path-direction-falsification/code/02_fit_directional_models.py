#!/usr/bin/env python3
"""Compare forward and reverse-compatible RT-ATAC-RNA path orderings.

The models assess asymmetry and falsification diagnostics; they do not identify
the true causal DAG because assays are unpaired and CPP8 loss can act in parallel.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import dmatrices
from scipy.stats import wilcoxon
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


SEED = 20260826
ORDERINGS = ("forward_rt_atac_rna", "reverse_atac_rt_rna", "reverse_rna_atac_rt")
FORMULAS = {
    "f_a": "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + z_log_n_ocr + C(chr)",
    "f_y": "z_rna_logfc ~ z_delta_atac + z_delta_rt + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
    "r1_a": "z_delta_rt ~ z_delta_atac + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
    "r2_a": "z_delta_atac ~ z_rna_logfc + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
    "r2_x": "z_delta_rt ~ z_delta_atac + z_rna_logfc + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--stage-contrasts", required=True)
    p.add_argument("--scale-contrasts", required=True)
    p.add_argument("--offsets", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=2000)
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def bootstrap_p(values: np.ndarray) -> float:
    x = values[np.isfinite(values)]
    lower = (np.sum(x <= 0) + 1) / (len(x) + 1)
    upper = (np.sum(x >= 0) + 1) / (len(x) + 1)
    return float(min(1.0, 2 * min(lower, upper)))


def sufficient(formula: str, data: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    aligned = data.loc[x.index, "block"].astype(str).to_numpy()
    lookup = {b: i for i, b in enumerate(blocks)}
    xa, ya = x.to_numpy(float), y.iloc[:, 0].to_numpy(float)
    xtx = np.zeros((len(blocks), x.shape[1], x.shape[1]))
    xty = np.zeros((len(blocks), x.shape[1]))
    for block in np.unique(aligned):
        take = aligned == block
        bi = lookup[block]
        xtx[bi] = xa[take].T @ xa[take]
        xty[bi] = xa[take].T @ ya[take]
    return {"xtx": xtx, "xty": xty, "columns": list(x.columns)}


def solve(suff: dict[str, object], sampled: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(suff["xtx"][sampled].sum(axis=0)) @ suff["xty"][sampled].sum(axis=0)


def coef(suff: dict[str, object], beta: np.ndarray, name: str) -> float:
    return float(beta[suff["columns"].index(name)])


def prepare(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["log_n_ocr"] = np.log1p(d["n_promoter_ocr"])
    needed = ["delta_rt", "delta_atac", "rna_logfc", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]
    d = d.dropna(subset=needed).copy()
    for col in needed:
        d[f"z_{col}"] = zscore(d[col])
    return d


def path_values(suffs: dict[str, dict[str, object]], sampled: np.ndarray) -> dict[str, dict[str, float]]:
    betas = {key: solve(suff, sampled) for key, suff in suffs.items()}
    f_a = coef(suffs["f_a"], betas["f_a"], "z_delta_rt")
    f_b = coef(suffs["f_y"], betas["f_y"], "z_delta_atac")
    r1_a = coef(suffs["r1_a"], betas["r1_a"], "z_delta_atac")
    r1_b = coef(suffs["f_y"], betas["f_y"], "z_delta_rt")
    r2_a = coef(suffs["r2_a"], betas["r2_a"], "z_rna_logfc")
    r2_b = coef(suffs["r2_x"], betas["r2_x"], "z_delta_atac")
    return {
        "forward_rt_atac_rna": {"a": f_a, "b": f_b, "indirect": f_a * f_b},
        "reverse_atac_rt_rna": {"a": r1_a, "b": r1_b, "indirect": r1_a * r1_b},
        "reverse_rna_atac_rt": {"a": r2_a, "b": r2_b, "indirect": r2_a * r2_b},
    }


RNA_MODELS = {
    "baseline": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"],
    "baseline_plus_rt": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_rt"],
    "baseline_plus_atac": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_atac"],
    "baseline_plus_rt_atac": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_rt", "delta_atac"],
}
RT_MODELS = {
    "baseline": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"],
    "baseline_plus_atac": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_atac"],
    "baseline_plus_rna": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "rna_logfc"],
    "baseline_plus_atac_rna": ["wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr", "delta_atac", "rna_logfc"],
}


def loco(d: pd.DataFrame, outcome: str, models: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    folds = []
    predictions = []
    for model_name, predictors in models.items():
        pred = np.full(len(d), np.nan)
        for chrom in sorted(d["chr"].unique()):
            test = d["chr"] == chrom
            model = make_pipeline(StandardScaler(), LinearRegression())
            model.fit(d.loc[~test, predictors], d.loc[~test, outcome])
            pred[test] = model.predict(d.loc[test, predictors])
            folds.append({"model": model_name, "left_out_chromosome": chrom, "mse": mean_squared_error(d.loc[test, outcome], pred[test]), "r2": r2_score(d.loc[test, outcome], pred[test]), "n_test": int(test.sum())})
        predictions.append(pd.DataFrame({"model": model_name, "observed": d[outcome], "predicted": pred}))
    fold_df = pd.DataFrame(folds)
    pred_df = pd.concat(predictions, ignore_index=True)
    summary = []
    for model_name in models:
        q = pred_df[pred_df["model"] == model_name]
        summary.append({"model": model_name, "pooled_loco_r2": r2_score(q["observed"], q["predicted"]), "pooled_loco_rmse": float(np.sqrt(mean_squared_error(q["observed"], q["predicted"])))})
    summary = pd.DataFrame(summary)
    baseline = float(summary.loc[summary["model"] == "baseline", "pooled_loco_r2"].iloc[0])
    summary["delta_r2_vs_baseline"] = summary["pooled_loco_r2"] - baseline
    return summary, fold_df


def paired_comparison(folds: pd.DataFrame, model_a: str, model_b: str, label: str) -> dict[str, object]:
    a = folds[folds["model"] == model_a].set_index("left_out_chromosome")["mse"]
    b = folds[folds["model"] == model_b].set_index("left_out_chromosome")["mse"]
    diff = a - b
    try:
        p = float(wilcoxon(diff, alternative="greater").pvalue)
    except ValueError:
        p = np.nan
    return {"comparison": label, "mean_mse_reduction": float(diff.mean()), "chromosome_wilcoxon_p": p}


def circular_shift_rna(d: pd.DataFrame) -> pd.Series:
    shifted = pd.Series(index=d.index, dtype=float)
    for _, idx in d.groupby("chr", sort=False).groups.items():
        ordered = np.asarray(sorted(idx, key=lambda i: d.at[i, "tss"]), dtype=int)
        shift = max(1, len(ordered) // 5)
        shifted.loc[ordered] = np.roll(d.loc[ordered, "rna_logfc"].to_numpy(float), shift)
    return shifted


def fit_forward_indirect(d: pd.DataFrame, nboot: int, rng: np.random.Generator) -> tuple[float, float, float, float]:
    d = prepare(d)
    blocks = np.array(sorted(d["block"].astype(str).unique()))
    sa = sufficient(FORMULAS["f_a"], d, blocks)
    sy = sufficient(FORMULAS["f_y"], d, blocks)
    ba = solve(sa, np.arange(len(blocks)))
    by = solve(sy, np.arange(len(blocks)))
    point = coef(sa, ba, "z_delta_rt") * coef(sy, by, "z_delta_atac")
    boots = np.empty(nboot)
    for i in range(nboot):
        sampled = rng.integers(0, len(blocks), len(blocks))
        ba, by = solve(sa, sampled), solve(sy, sampled)
        boots[i] = coef(sa, ba, "z_delta_rt") * coef(sy, by, "z_delta_atac")
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi), bootstrap_p(boots)


def main() -> None:
    args = parse_args()
    outdir, plotdir = Path(args.output_dir).resolve(), Path(args.plot_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if (task_root not in outdir.parents and outdir != task_root) or (task_root not in plotdir.parents and plotdir != task_root):
        raise ValueError("Refusing output outside task root")
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.gene_metrics, sep="\t")

    direction_rows, contrast_rows, loco_rows, loco_compare_rows = [], [], [], []
    group_data = {}
    for keys, raw in metrics.groupby(["allele", "stage", "promoter_window_bp"], sort=False):
        allele, stage, window = keys
        d = prepare(raw.reset_index(drop=True))
        group_data[keys] = d
        blocks = np.array(sorted(d["block"].astype(str).unique()))
        suffs = {key: sufficient(formula, d, blocks) for key, formula in FORMULAS.items()}
        point = path_values(suffs, np.arange(len(blocks)))
        boots = {ordering: {part: np.empty(args.bootstrap) for part in ("a", "b", "indirect")} for ordering in ORDERINGS}
        for i in range(args.bootstrap):
            vals = path_values(suffs, rng.integers(0, len(blocks), len(blocks)))
            for ordering in ORDERINGS:
                for part in ("a", "b", "indirect"):
                    boots[ordering][part][i] = vals[ordering][part]
        for ordering in ORDERINGS:
            row = {"allele": allele, "stage": stage, "promoter_window_bp": int(window), "ordering": ordering, "n_genes": len(d), "n_blocks": len(blocks)}
            for part in ("a", "b", "indirect"):
                lo, hi = np.percentile(boots[ordering][part], [2.5, 97.5])
                row[f"path_{part}"] = point[ordering][part]
                row[f"{part}_boot_ci_low"] = float(lo)
                row[f"{part}_boot_ci_high"] = float(hi)
                row[f"{part}_bootstrap_p"] = bootstrap_p(boots[ordering][part])
            direction_rows.append(row)
        for reverse in ("reverse_atac_rt_rna", "reverse_rna_atac_rt"):
            diff = boots["forward_rt_atac_rna"]["indirect"] - boots[reverse]["indirect"]
            lo, hi = np.percentile(diff, [2.5, 97.5])
            contrast_rows.append({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "contrast": f"forward_minus_{reverse}", "estimate": point["forward_rt_atac_rna"]["indirect"] - point[reverse]["indirect"], "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(diff)})

        for outcome_name, outcome, models in [("rna", "rna_logfc", RNA_MODELS), ("rt", "delta_rt", RT_MODELS)]:
            summary, folds = loco(d, outcome, models)
            for _, row in summary.iterrows():
                loco_rows.append({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "outcome": outcome_name, **row.to_dict()})
            if outcome_name == "rna":
                comps = [
                    paired_comparison(folds, "baseline_plus_rt", "baseline_plus_rt_atac", "add_atac_given_rt"),
                    paired_comparison(folds, "baseline_plus_atac", "baseline_plus_rt_atac", "add_rt_given_atac"),
                ]
            else:
                comps = [
                    paired_comparison(folds, "baseline_plus_rna", "baseline_plus_atac_rna", "add_atac_given_rna"),
                    paired_comparison(folds, "baseline_plus_atac", "baseline_plus_atac_rna", "add_rna_given_atac"),
                ]
            for row in comps:
                loco_compare_rows.append({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "outcome": outcome_name, **row})

    directions = pd.DataFrame(direction_rows)
    directions["indirect_fdr_across_36_ordering_tests"] = multipletests(directions["indirect_bootstrap_p"], method="fdr_bh")[1]
    directions.to_csv(outdir / "all_direction_path_effects_013.tsv", sep="\t", index=False)
    direction_contrasts = pd.DataFrame(contrast_rows)
    direction_contrasts["fdr_within_contrast"] = direction_contrasts.groupby("contrast")["bootstrap_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    direction_contrasts.to_csv(outdir / "forward_vs_reverse_contrasts_013.tsv", sep="\t", index=False)
    loco_df = pd.DataFrame(loco_rows)
    loco_comp = pd.DataFrame(loco_compare_rows)
    loco_df.to_csv(outdir / "directional_loco_summary_013.tsv", sep="\t", index=False)
    loco_comp.to_csv(outdir / "directional_loco_comparisons_013.tsv", sep="\t", index=False)

    negative_rows = []
    mismatch = {"ES": "LS", "MS": "ES", "LS": "ES"}
    for (allele, stage, window), d in group_data.items():
        other = group_data[(allele, mismatch[stage], window)][["gene_id", "delta_atac", "wt_log_atac"]].rename(columns={"delta_atac": "mismatch_delta_atac", "wt_log_atac": "mismatch_wt_log_atac"})
        wrong = d.drop(columns=["delta_atac", "wt_log_atac"]).merge(other, on="gene_id").rename(columns={"mismatch_delta_atac": "delta_atac", "mismatch_wt_log_atac": "wt_log_atac"})
        point, lo, hi, p = fit_forward_indirect(wrong, min(args.bootstrap, 1000), rng)
        negative_rows.append({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "negative_control": f"wrong_stage_{mismatch[stage]}", "indirect_effect": point, "bootstrap_ci_low": lo, "bootstrap_ci_high": hi, "bootstrap_p": p})
        shifted = d.copy()
        shifted["rna_logfc"] = circular_shift_rna(shifted)
        point, lo, hi, p = fit_forward_indirect(shifted, min(args.bootstrap, 1000), rng)
        negative_rows.append({"allele": allele, "stage": stage, "promoter_window_bp": int(window), "negative_control": "chr_circular_shift_rna", "indirect_effect": point, "bootstrap_ci_low": lo, "bootstrap_ci_high": hi, "bootstrap_p": p})
    negatives = pd.DataFrame(negative_rows)
    negatives["fdr_within_negative_control"] = negatives.groupby("negative_control")["bootstrap_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    negatives.to_csv(outdir / "negative_control_paths_013.tsv", sep="\t", index=False)

    stage = pd.read_csv(args.stage_contrasts, sep="\t")
    scale = pd.read_csv(args.scale_contrasts, sep="\t")
    offsets = pd.read_csv(args.offsets, sep="\t")
    evidence_rows = []
    for allele in ("cpp8_1", "cpp8_3"):
        for window in (1000, 3000):
            f = directions[(directions["allele"] == allele) & (directions["stage"] == "MS") & (directions["promoter_window_bp"] == window) & (directions["ordering"] == "forward_rt_atac_rna")].iloc[0]
            reverse = directions[(directions["allele"] == allele) & (directions["stage"] == "MS") & (directions["promoter_window_bp"] == window) & (directions["ordering"] != "forward_rt_atac_rna")]
            evidence_rows.append({"diagnostic": "two_allele_forward_path", "allele": allele, "promoter_window_bp": window, "status": "supports_forward_compatibility" if f["indirect_boot_ci_low"] > 0 else "not_supported", "evidence": f"IF={f['path_indirect']:.4f}, 95% CI {f['indirect_boot_ci_low']:.4f} to {f['indirect_boot_ci_high']:.4f}"})
            evidence_rows.append({"diagnostic": "reverse_orderings", "allele": allele, "promoter_window_bp": window, "status": "reverse_warning" if (reverse["indirect_boot_ci_low"] > 0).any() | (reverse["indirect_boot_ci_high"] < 0).any() else "reverse_not_stable", "evidence": "; ".join(f"{r.ordering}={r.path_indirect:.4f}" for r in reverse.itertuples())})
        st = stage[(stage["allele"] == allele) & (stage["contrast"] == "global_ms_dominance")].iloc[0]
        evidence_rows.append({"diagnostic": "temporal_structure", "allele": allele, "promoter_window_bp": 0, "status": "global_ms_supported" if st["bootstrap_ci_low"] > 0 else "not_supported", "evidence": f"MS dominance={st['estimate']:.4f}"})
        sc = scale[(scale["allele"] == allele) & (scale["stage"] == "MS") & (scale["contrast"] == "local_dominance")].iloc[0]
        off = offsets[(offsets["allele"] == allele) & (offsets["stage"] == "MS")]
        peak = int(off.loc[off["standardized_beta_delta_rt"].idxmax(), "offset_kb"])
        evidence_rows.append({"diagnostic": "spatial_structure", "allele": allele, "promoter_window_bp": 0, "status": "local_peak" if peak == 0 else "offset_peak", "evidence": f"local dominance={sc['estimate']:.4f}; offset peak={peak} kb"})
    pd.DataFrame(evidence_rows).to_csv(outdir / "directional_evidence_matrix_013.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    primary = directions[(directions["stage"] == "MS")].copy()
    primary["label"] = primary["allele"].str.replace("cpp8_", "oscpp8-") + " ±" + (primary["promoter_window_bp"] // 1000).astype(str) + "kb"
    labels = list(dict.fromkeys(primary["label"]))
    ypos = {label: i for i, label in enumerate(labels)}
    colors = {"forward_rt_atac_rna": "#228833", "reverse_atac_rt_rna": "#AA4499", "reverse_rna_atac_rt": "#4477AA"}
    offsets_y = {"forward_rt_atac_rna": -0.18, "reverse_atac_rt_rna": 0.0, "reverse_rna_atac_rt": 0.18}
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for ordering, q in primary.groupby("ordering"):
        y = [ypos[x] + offsets_y[ordering] for x in q["label"]]
        ax.errorbar(q["path_indirect"], y, xerr=[q["path_indirect"] - q["indirect_boot_ci_low"], q["indirect_boot_ci_high"] - q["path_indirect"]], fmt="o", capsize=3, color=colors[ordering], label=ordering.replace("_", " "))
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("Standardized product path")
    ax.set_title("Forward and reverse-compatible MS paths")
    ax.legend(frameon=False, fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "forward_reverse_path_forest_013.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "forward_reverse_path_forest_013.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Directional paths")
    print(directions.to_string(index=False))
    print("\nForward-minus-reverse contrasts")
    print(direction_contrasts.to_string(index=False))
    print("\nNegative controls")
    print(negatives.to_string(index=False))


if __name__ == "__main__":
    main()
