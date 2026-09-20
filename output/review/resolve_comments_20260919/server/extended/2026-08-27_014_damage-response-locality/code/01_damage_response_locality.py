#!/usr/bin/env python3
"""Test whether stress/cell-cycle transcriptional programs localize to CPP8 RT/ATAC disruption.

The primary universe contains genes detected in the joint RNA analysis with at
least one OCR in the TSS +/-1-kb promoter.  The primary gene sets are fixed
before testing: the manuscript DNA-repair panel, mapped rice orthologues of
Arabidopsis SOG1 targets, and a curated rice cell-cycle set.

Inference uses 1-Mb genomic-block bootstrap.  A covariate-stratified label
permutation is reported as a secondary matched null.  This is a locality test,
not individual-level causal mediation.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
from patsy import dmatrices
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


SEED = 20260826
ALLELES = ("cpp8_1", "cpp8_3")
STAGES = ("ES", "MS", "LS")
GENE_SET_ORDER = ("manuscript_dna_repair", "sog1_targets", "cell_cycle")
GENE_SET_LABELS = {
    "manuscript_dna_repair": "DNA-repair panel",
    "sog1_targets": "SOG1 targets",
    "cell_cycle": "Cell-cycle genes",
}
COLORS = {
    "manuscript_dna_repair": "#9E3D45",
    "sog1_targets": "#D17B39",
    "cell_cycle": "#3C6E9E",
}
HOTSPOT_THRESHOLD = 0.2632919444
HOTSPOT_TOP10_THRESHOLD = 0.2046691


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--gene-tss", required=True)
    p.add_argument("--dna-repair", required=True)
    p.add_argument("--cell-cycle", required=True)
    p.add_argument("--sog1", required=True)
    p.add_argument("--rice2ara", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--permutations", type=int, default=5000)
    return p.parse_args()


def normalize_gene_id(value: object) -> str:
    x = str(value).strip().replace("\ufeff", "")
    return re.sub(r"^(LOC_Os\d{2})G", r"\1g", x)


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=series.index)
    return (x - x.mean()) / sd


def scaled_magnitude(series: pd.Series) -> pd.Series:
    """Log-transform a non-negative magnitude and standardize it."""
    x = pd.to_numeric(series, errors="coerce").clip(lower=0)
    positive = x[x > 0]
    scale = float(positive.median()) if len(positive) else 1.0
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    return zscore(np.log1p(x / scale))


def conservative_consensus(a: pd.Series, b: pd.Series) -> pd.Series:
    av = pd.to_numeric(a, errors="coerce").to_numpy(float)
    bv = pd.to_numeric(b, errors="coerce").to_numpy(float)
    same = np.isfinite(av) & np.isfinite(bv) & ((av * bv) > 0)
    out = np.zeros(len(av), dtype=float)
    out[same] = np.minimum(np.abs(av[same]), np.abs(bv[same]))
    out[~(np.isfinite(av) & np.isfinite(bv))] = np.nan
    return pd.Series(out, index=a.index)


def read_gene_sets(args: argparse.Namespace) -> tuple[dict[str, set[str]], pd.DataFrame]:
    dna = pd.read_csv(args.dna_repair)
    dna_set = {normalize_gene_id(x) for x in dna.iloc[:, 0].dropna()}

    cell = pd.read_csv(args.cell_cycle)
    if "msu" not in cell.columns:
        raise ValueError("Cell-cycle table lacks the msu column")
    cell_set = {normalize_gene_id(x) for x in cell["msu"].dropna()}

    sog = pd.read_csv(args.sog1, header=None)
    at_targets = {str(x).strip() for x in sog.iloc[:, 0].dropna()}
    trans = pd.read_csv(args.rice2ara)
    trans["rice"] = trans["rice"].map(normalize_gene_id)
    trans["tair"] = trans["tair"].astype(str).str.strip()
    sog_set = set(trans.loc[trans["tair"].isin(at_targets), "rice"].dropna())

    sets = {
        "manuscript_dna_repair": dna_set,
        "sog1_targets": sog_set,
        "cell_cycle": cell_set,
    }
    membership = sorted(set().union(*sets.values()))
    rows = []
    for gene in membership:
        rows.append({"gene_id": gene, **{name: gene in genes for name, genes in sets.items()}})
    return sets, pd.DataFrame(rows)


def build_window_table(metrics: pd.DataFrame, gene_tss: pd.DataFrame, window: int) -> pd.DataFrame:
    d = metrics[metrics["promoter_window_bp"] == window].copy()
    key = ["gene_id", "chr", "tss", "strand", "block"]
    base = d[(d["allele"] == "cpp8_1") & (d["stage"] == "MS")][
        key + ["n_promoter_ocr", "wt_wrt", "baseline_rna"]
    ].drop_duplicates("gene_id")

    for stage in STAGES:
        wt = d[(d["allele"] == "cpp8_1") & (d["stage"] == stage)][
            ["gene_id", "wt_log_atac"]
        ].drop_duplicates("gene_id").rename(columns={"wt_log_atac": f"wt_log_atac_{stage}"})
        base = base.merge(wt, on="gene_id", how="left", validate="one_to_one")

    for allele in ALLELES:
        a0 = d[(d["allele"] == allele) & (d["stage"] == "MS")][
            ["gene_id", "delta_rt", "rna_logfc"]
        ].drop_duplicates("gene_id").rename(
            columns={"delta_rt": f"delta_rt_{allele}", "rna_logfc": f"rna_logfc_{allele}"}
        )
        base = base.merge(a0, on="gene_id", how="left", validate="one_to_one")
        for stage in STAGES:
            a = d[(d["allele"] == allele) & (d["stage"] == stage)][
                ["gene_id", "delta_atac"]
            ].drop_duplicates("gene_id").rename(columns={"delta_atac": f"delta_atac_{allele}_{stage}"})
            base = base.merge(a, on="gene_id", how="left", validate="one_to_one")

    lengths = gene_tss[["gene_id", "start", "end"]].drop_duplicates("gene_id").copy()
    lengths["gene_length"] = pd.to_numeric(lengths["end"], errors="coerce") - pd.to_numeric(
        lengths["start"], errors="coerce"
    ) + 1
    base = base.merge(lengths[["gene_id", "gene_length"]], on="gene_id", how="left", validate="one_to_one")

    base["rt_consensus_mag"] = conservative_consensus(
        base["delta_rt_cpp8_1"], base["delta_rt_cpp8_3"]
    )
    base["rt_hotspot_top5"] = (base["rt_consensus_mag"] >= HOTSPOT_THRESHOLD).astype(float)
    base["rt_hotspot_top10"] = (base["rt_consensus_mag"] >= HOTSPOT_TOP10_THRESHOLD).astype(float)
    base["rt_gene_relative_top5"] = (
        base["rt_consensus_mag"] >= base["rt_consensus_mag"].quantile(0.95)
    ).astype(float)
    base["rna_consensus_mag"] = conservative_consensus(
        base["rna_logfc_cpp8_1"], base["rna_logfc_cpp8_3"]
    )
    for allele in ALLELES:
        base[f"rt_abs_{allele}"] = base[f"delta_rt_{allele}"].abs()
        base[f"rna_abs_{allele}"] = base[f"rna_logfc_{allele}"].abs()
    for stage in STAGES:
        base[f"atac_consensus_mag_{stage}"] = conservative_consensus(
            base[f"delta_atac_cpp8_1_{stage}"], base[f"delta_atac_cpp8_3_{stage}"]
        )
        for allele in ALLELES:
            base[f"atac_abs_{allele}_{stage}"] = base[f"delta_atac_{allele}_{stage}"].abs()

    base["log_gene_length"] = np.log1p(pd.to_numeric(base["gene_length"], errors="coerce"))
    base["log_n_ocr"] = np.log1p(pd.to_numeric(base["n_promoter_ocr"], errors="coerce"))
    base["z_wt_wrt"] = zscore(base["wt_wrt"])
    base["z_baseline_rna"] = zscore(base["baseline_rna"])
    base["z_log_gene_length"] = zscore(base["log_gene_length"])
    base["z_log_n_ocr"] = zscore(base["log_n_ocr"])
    base["z_wt_log_atac"] = zscore(base["wt_log_atac_MS"])
    for col in [
        "rt_consensus_mag",
        "rna_consensus_mag",
        *[f"rt_abs_{a}" for a in ALLELES],
        *[f"rna_abs_{a}" for a in ALLELES],
        *[f"atac_consensus_mag_{s}" for s in STAGES],
        *[f"atac_abs_{a}_{s}" for a in ALLELES for s in STAGES],
    ]:
        base[f"z_{col}"] = scaled_magnitude(base[col])

    # Strata used only for the secondary covariate-matched label permutation.
    base["match_wrt"] = pd.qcut(base["wt_wrt"], 5, labels=False, duplicates="drop")
    base["match_rna"] = pd.qcut(base["baseline_rna"], 5, labels=False, duplicates="drop")
    base["match_stratum"] = (
        base["chr"].astype(str)
        + ":"
        + base["match_wrt"].astype("Int64").astype(str)
        + ":"
        + base["match_rna"].astype("Int64").astype(str)
    )
    return base


def sufficient_stats(formula: str, data: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    block_series = data.loc[x.index, "block"].astype(str)
    lookup = {b: i for i, b in enumerate(blocks)}
    xa = x.to_numpy(float)
    ya = y.iloc[:, 0].to_numpy(float)
    xtx = np.zeros((len(blocks), xa.shape[1], xa.shape[1]), dtype=float)
    xty = np.zeros((len(blocks), xa.shape[1]), dtype=float)
    for block, idx in block_series.groupby(block_series).groups.items():
        positions = x.index.get_indexer(idx)
        bi = lookup[block]
        xb = xa[positions]
        yb = ya[positions]
        xtx[bi] = xb.T @ xb
        xty[bi] = xb.T @ yb
    return {"xtx": xtx, "xty": xty, "columns": list(x.columns), "n": len(x)}


def solve_sufficient(suff: dict[str, object], sampled: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(suff["xtx"][sampled].sum(axis=0)) @ suff["xty"][sampled].sum(axis=0)


def two_sided_bootstrap_p(values: np.ndarray) -> float:
    x = values[np.isfinite(values)]
    if not len(x):
        return np.nan
    lower = (np.sum(x <= 0) + 1) / (len(x) + 1)
    upper = (np.sum(x >= 0) + 1) / (len(x) + 1)
    return float(min(1.0, 2 * min(lower, upper)))


def matched_permutation_p(
    data: pd.DataFrame,
    outcome: str,
    set_col: str,
    nperm: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    cov_formula = (
        f"{outcome} ~ z_wt_wrt + z_baseline_rna + z_log_gene_length + "
        "z_log_n_ocr + z_wt_log_atac + C(chr)"
    )
    fit = smf.ols(cov_formula, data=data).fit()
    q = data.loc[fit.model.data.row_labels].copy()
    q["residual"] = fit.resid
    target = q[set_col].astype(bool).to_numpy()
    point = float(q.loc[target, "residual"].mean() - q.loc[~target, "residual"].mean())
    strata = []
    for _, idx in q.groupby("match_stratum", observed=True).groups.items():
        pos = q.index.get_indexer(idx)
        k = int(target[pos].sum())
        if k:
            strata.append((pos, k))
    residual = q["residual"].to_numpy(float)
    stats = np.empty(nperm, dtype=float)
    all_idx = np.arange(len(q))
    for i in range(nperm):
        selected_parts = []
        for pos, k in strata:
            selected_parts.append(rng.choice(pos, size=k, replace=False))
        selected = np.concatenate(selected_parts) if selected_parts else np.array([], dtype=int)
        mask = np.zeros(len(q), dtype=bool)
        mask[selected] = True
        stats[i] = residual[mask].mean() - residual[~mask].mean()
    p = (np.sum(np.abs(stats) >= abs(point)) + 1) / (nperm + 1)
    return point, float(p)


def fit_endpoint(
    data: pd.DataFrame,
    outcome: str,
    set_name: str,
    endpoint_label: str,
    endpoint_type: str,
    window: int,
    family: str,
    nboot: int,
    nperm: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    set_col = f"set_{set_name}"
    needed = [
        outcome,
        set_col,
        "block",
        "chr",
        "z_wt_wrt",
        "z_baseline_rna",
        "z_log_gene_length",
        "z_log_n_ocr",
        "z_wt_log_atac",
        "match_stratum",
    ]
    q = data.dropna(subset=needed).copy().reset_index(drop=True)
    coefficient = set_col
    formula = (
        f"{outcome} ~ {set_col} + z_wt_wrt + z_baseline_rna + z_log_gene_length + "
        "z_log_n_ocr + z_wt_log_atac + C(chr)"
    )
    cluster_fit = smf.ols(formula, data=q).fit(cov_type="cluster", cov_kwds={"groups": q["block"]})
    blocks = np.sort(q["block"].astype(str).unique())
    suff = sufficient_stats(formula, q, blocks)
    coef_idx = suff["columns"].index(coefficient)
    point = float(cluster_fit.params[coefficient])
    boots = np.empty(nboot, dtype=float)
    for i in range(nboot):
        sampled = rng.integers(0, len(blocks), len(blocks))
        boots[i] = solve_sufficient(suff, sampled)[coef_idx]
    ci_low, ci_high = np.nanpercentile(boots, [2.5, 97.5])
    matched_diff, perm_p = matched_permutation_p(q, outcome, set_col, nperm, rng)

    target = q[set_col].astype(bool)
    raw_col = endpoint_label
    if raw_col not in q.columns:
        raw_col = outcome
    return {
        "family": family,
        "window_bp": window,
        "gene_set": set_name,
        "gene_set_label": GENE_SET_LABELS[set_name],
        "endpoint": endpoint_label,
        "endpoint_type": endpoint_type,
        "outcome_column": outcome,
        "n_universe": len(q),
        "n_set": int(target.sum()),
        "n_blocks": len(blocks),
        "adjusted_effect": point,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "cluster_p": float(cluster_fit.pvalues[coefficient]),
        "bootstrap_p": two_sided_bootstrap_p(boots),
        "matched_residual_difference": matched_diff,
        "matched_permutation_p": perm_p,
        "set_raw_mean": float(q.loc[target, raw_col].mean()),
        "control_raw_mean": float(q.loc[~target, raw_col].mean()),
        "set_raw_median": float(q.loc[target, raw_col].median()),
        "control_raw_median": float(q.loc[~target, raw_col].median()),
    }


def run_models(
    tables: dict[int, pd.DataFrame],
    nboot: int,
    nperm: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    primary = [
        ("z_rt_consensus_mag", "rt_consensus_mag", "standardized_magnitude"),
        ("z_atac_consensus_mag_MS", "atac_consensus_mag_MS", "standardized_magnitude"),
        ("rt_hotspot_top5", "rt_hotspot_top5", "risk_difference"),
    ]
    for set_name in GENE_SET_ORDER:
        for outcome, label, typ in primary:
            rows.append(
                fit_endpoint(
                    tables[1000], outcome, set_name, label, typ, 1000, "primary", nboot, nperm, rng
                )
            )

    # Predefined robustness: wider promoter universe, ES/LS accessibility,
    # and allele-specific magnitude estimates.
    for set_name in GENE_SET_ORDER:
        for outcome, label, typ in primary:
            rows.append(
                fit_endpoint(
                    tables[3000], outcome, set_name, label, typ, 3000, "window_sensitivity", nboot, nperm, rng
                )
            )
        for outcome, label in [
            ("rt_hotspot_top10", "rt_hotspot_top10"),
            ("rt_gene_relative_top5", "rt_gene_relative_top5"),
        ]:
            rows.append(
                fit_endpoint(
                    tables[1000],
                    outcome,
                    set_name,
                    label,
                    "risk_difference",
                    1000,
                    "hotspot_sensitivity",
                    nboot,
                    nperm,
                    rng,
                )
            )
        for stage in ("ES", "LS"):
            rows.append(
                fit_endpoint(
                    tables[1000],
                    f"z_atac_consensus_mag_{stage}",
                    set_name,
                    f"atac_consensus_mag_{stage}",
                    "standardized_magnitude",
                    1000,
                    "stage_sensitivity",
                    nboot,
                    nperm,
                    rng,
                )
            )
        for allele in ALLELES:
            rows.append(
                fit_endpoint(
                    tables[1000],
                    f"z_rt_abs_{allele}",
                    set_name,
                    f"rt_abs_{allele}",
                    "standardized_magnitude",
                    1000,
                    "allele_sensitivity",
                    nboot,
                    nperm,
                    rng,
                )
            )
            for stage in STAGES:
                rows.append(
                    fit_endpoint(
                        tables[1000],
                        f"z_atac_abs_{allele}_{stage}",
                        set_name,
                        f"atac_abs_{allele}_{stage}",
                        "standardized_magnitude",
                        1000,
                        "allele_stage_sensitivity",
                        nboot,
                        nperm,
                        rng,
                    )
                )
        # RNA is a descriptive positive control.  It is outcome-selected for
        # the manuscript DNA-repair panel and is never part of the primary FDR.
        rows.append(
            fit_endpoint(
                tables[1000],
                "z_rna_consensus_mag",
                set_name,
                "rna_consensus_mag",
                "standardized_magnitude",
                1000,
                "rna_descriptive",
                nboot,
                nperm,
                rng,
            )
        )
    result = pd.DataFrame(rows)
    result["bootstrap_fdr"] = np.nan
    result["permutation_fdr"] = np.nan
    result["inference_p_source"] = np.where(
        result["endpoint_type"] == "risk_difference", "matched_permutation", "block_bootstrap"
    )
    result["inference_p"] = np.where(
        result["endpoint_type"] == "risk_difference",
        result["matched_permutation_p"],
        result["bootstrap_p"],
    )
    result["inference_fdr"] = np.nan
    for family, idx in result.groupby("family").groups.items():
        result.loc[idx, "bootstrap_fdr"] = multipletests(
            result.loc[idx, "bootstrap_p"].fillna(1.0), method="fdr_bh"
        )[1]
        result.loc[idx, "permutation_fdr"] = multipletests(
            result.loc[idx, "matched_permutation_p"].fillna(1.0), method="fdr_bh"
        )[1]
        result.loc[idx, "inference_fdr"] = multipletests(
            result.loc[idx, "inference_p"].fillna(1.0), method="fdr_bh"
        )[1]
    return result


def plot_results(results: pd.DataFrame, plot_dir: Path) -> None:
    primary = results[(results["family"] == "primary") & (results["window_bp"] == 1000)].copy()
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2), gridspec_kw={"width_ratios": [1, 1, 1]})
    endpoints = [
        ("rt_consensus_mag", "Two-allele consensus RT\ninstability", "Adjusted effect (s.d.)"),
        ("atac_consensus_mag_MS", "Two-allele consensus MS\naccessibility disruption", "Adjusted effect (s.d.)"),
        ("rt_hotspot_top5", "Top-5% RT-instability\nhotspot occupancy", "Fraction of genes"),
    ]
    y_map = {name: i for i, name in enumerate(reversed(GENE_SET_ORDER))}
    for ax, (endpoint, title, xlabel) in zip(axes, endpoints):
        q = primary[primary["endpoint"] == endpoint]
        if endpoint != "rt_hotspot_top5":
            ax.axvline(0, color="#6F6F6F", lw=0.8, ls="--")
        for _, row in q.iterrows():
            y = y_map[row["gene_set"]]
            if endpoint == "rt_hotspot_top5":
                ax.plot(
                    [row["set_raw_mean"], row["control_raw_mean"]], [y, y],
                    color="#B8B8B8", lw=1.3, zorder=1,
                )
                ax.scatter(row["control_raw_mean"], y, marker="D", s=33, facecolor="white", edgecolor="#666666", zorder=2)
                ax.scatter(row["set_raw_mean"], y, s=45, color=COLORS[row["gene_set"]], zorder=3)
                ax.text(
                    max(row["set_raw_mean"], row["control_raw_mean"]) + 0.00025,
                    y,
                    f"matched P={row['matched_permutation_p']:.2f}",
                    va="center",
                    fontsize=7.5,
                    color="#555555",
                )
            else:
                ax.errorbar(
                    row["adjusted_effect"],
                    y,
                    xerr=[[row["adjusted_effect"] - row["ci_low"]], [row["ci_high"] - row["adjusted_effect"]]],
                    fmt="o",
                    ms=6,
                    color=COLORS[row["gene_set"]],
                    ecolor=COLORS[row["gene_set"]],
                    capsize=3,
                )
                if row["inference_fdr"] < 0.05:
                    ax.text(row["ci_high"], y + 0.12, "*", color=COLORS[row["gene_set"]], fontsize=11)
        ax.set_yticks(range(len(GENE_SET_ORDER)))
        ax.set_yticklabels([GENE_SET_LABELS[x] for x in reversed(GENE_SET_ORDER)] if ax is axes[0] else [])
        ax.set_title(title, fontsize=10.5, weight="bold")
        ax.set_xlabel(xlabel, fontsize=9.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8.5)
        if endpoint == "rt_hotspot_top5":
            ax.set_xlim(-0.00035, 0.0082)
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[x], markeredgecolor=COLORS[x], label=GENE_SET_LABELS[x])
        for x in GENE_SET_ORDER
    ]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Do DNA-damage and cell-cycle programmes localize to CPP8-dependent RT/ATAC disruption?", fontsize=12, weight="bold")
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    for ext in ("pdf", "svg", "png"):
        kwargs = {"dpi": 300} if ext == "png" else {}
        fig.savefig(plot_dir / f"damage_response_locality_014.{ext}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    stage = results[(results["family"].isin(["primary", "stage_sensitivity"])) & (results["window_bp"] == 1000)].copy()
    stage = stage[stage["endpoint"].isin([f"atac_consensus_mag_{s}" for s in STAGES])]
    matrix = stage.pivot(index="gene_set", columns="endpoint", values="adjusted_effect").reindex(GENE_SET_ORDER)
    matrix = matrix.reindex(columns=[f"atac_consensus_mag_{s}" for s in STAGES])
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    vmax = max(0.1, float(np.nanmax(np.abs(matrix.to_numpy()))))
    im = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(3), STAGES)
    ax.set_yticks(range(3), [GENE_SET_LABELS[x] for x in GENE_SET_ORDER])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iloc[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title(
        "Adjusted enrichment of two-allele accessibility disruption\n(no stage-specific test survived FDR correction)",
        fontsize=10.0,
        weight="bold",
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cbar.set_label("Adjusted effect (s.d.)", fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        kwargs = {"dpi": 300} if ext == "png" else {}
        fig.savefig(plot_dir / f"damage_response_stage_heatmap_014.{ext}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    metrics = pd.read_csv(args.gene_metrics, sep="\t")
    metrics["gene_id"] = metrics["gene_id"].map(normalize_gene_id)
    gene_tss = pd.read_csv(args.gene_tss, sep="\t")
    gene_tss["gene_id"] = gene_tss["gene_id"].map(normalize_gene_id)
    gene_sets, membership = read_gene_sets(args)

    tables: dict[int, pd.DataFrame] = {}
    for window in (1000, 3000):
        table = build_window_table(metrics, gene_tss, window)
        for set_name, genes in gene_sets.items():
            table[f"set_{set_name}"] = table["gene_id"].isin(genes).astype(int)
        tables[window] = table
        table.to_csv(outdir / f"gene_universe_{window}bp_014.tsv.gz", sep="\t", index=False, compression="gzip")

    universe = set(tables[1000]["gene_id"])
    membership["in_primary_universe"] = membership["gene_id"].isin(universe)
    membership.to_csv(outdir / "gene_set_membership_014.tsv", sep="\t", index=False)

    overlap_rows = []
    for a in GENE_SET_ORDER:
        for b in GENE_SET_ORDER:
            sa = gene_sets[a] & universe
            sb = gene_sets[b] & universe
            overlap_rows.append(
                {
                    "set_a": a,
                    "set_b": b,
                    "n_a_in_universe": len(sa),
                    "n_b_in_universe": len(sb),
                    "intersection": len(sa & sb),
                    "jaccard": len(sa & sb) / len(sa | sb) if (sa | sb) else np.nan,
                }
            )
    pd.DataFrame(overlap_rows).to_csv(outdir / "gene_set_overlap_014.tsv", sep="\t", index=False)

    results = run_models(tables, args.bootstrap, args.permutations, rng)
    results.to_csv(outdir / "damage_response_locality_models_014.tsv", sep="\t", index=False)
    results[results["family"] == "primary"].to_csv(
        outdir / "damage_response_locality_primary_014.tsv", sep="\t", index=False
    )
    results[results["family"] != "primary"].to_csv(
        outdir / "damage_response_locality_sensitivity_014.tsv", sep="\t", index=False
    )
    plot_results(results, plotdir)

    primary = results[results["family"] == "primary"].copy()
    summary = {
        "seed": SEED,
        "bootstrap": args.bootstrap,
        "permutations": args.permutations,
        "hotspot_threshold": HOTSPOT_THRESHOLD,
        "primary_universe_n": len(tables[1000]),
        "sensitivity_universe_n": len(tables[3000]),
        "gene_set_sizes_in_primary_universe": {
            name: int(tables[1000][f"set_{name}"].sum()) for name in GENE_SET_ORDER
        },
        "primary_fdr_lt_0_05": int((primary["inference_fdr"] < 0.05).sum()),
        "primary_results": primary[
            [
                "gene_set", "endpoint", "adjusted_effect", "ci_low", "ci_high",
                "inference_p_source", "inference_p", "inference_fdr",
            ]
        ].to_dict(orient="records"),
    }
    (outdir / "damage_response_locality_summary_014.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
