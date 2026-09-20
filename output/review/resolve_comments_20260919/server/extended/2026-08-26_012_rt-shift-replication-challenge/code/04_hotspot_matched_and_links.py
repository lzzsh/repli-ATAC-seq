#!/usr/bin/env python3
"""Matched hotspot enrichment and links to ATAC/RNA perturbation magnitudes."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


SEED = 20260826
PRIMARY_FEATURES = {
    "log_nearest_gene_length": "Nearest-gene length",
    "wt_expression_50kb_log2cpm": "WT expression within 50 kb",
    "class_i_te_fraction": "Class I TE fraction",
    "class_ii_te_fraction": "Class II TE fraction",
    "gc_fraction": "GC fraction",
    "log_distance_to_centromere_proxy": "Distance to centromere-repeat proxy",
    "log_gene_count_50kb": "Gene density within 50 kb",
    "log_ocr_count_50kb": "OCR density within 50 kb",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True)
    p.add_argument("--ocr-metrics", required=True)
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=1000)
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


def quantile_label(s: pd.Series, q: int) -> pd.Series:
    return pd.qcut(s.rank(method="first"), q, labels=False, duplicates="drop")


def prepare_features(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["log_nearest_gene_length"] = np.log1p(d["nearest_gene_length"])
    d["log_distance_to_centromere_proxy"] = np.log1p(d["distance_to_centromere_proxy"])
    d["log_gene_count_50kb"] = np.log1p(d["gene_count_50kb"])
    d["log_ocr_count_50kb"] = np.log1p(d["ocr_count_50kb"])
    for col in PRIMARY_FEATURES:
        d[f"z_{col}"] = zscore(d[col])
    d["wrt_decile"] = quantile_label(d["WT_wrt"], 10)
    d["gc_quintile"] = quantile_label(d["gc_fraction"], 5)
    d["g1_quintile"] = quantile_label(d["log_wt_g1_coverage"], 5)
    gene = d["gene_body_fraction"] > 0
    te1 = d["class_i_te_fraction"] > 0
    te2 = d["class_ii_te_fraction"] > 0
    d["compartment"] = np.select([gene, te1 & te2, te1, te2], ["gene", "both_te", "class_i", "class_ii"], default="other")
    d["stratum"] = d["chr"].astype(str) + ":" + d["wrt_decile"].astype(str) + ":" + d["gc_quintile"].astype(str) + ":" + d["g1_quintile"].astype(str) + ":" + d["compartment"]
    return d


def matched_analysis(d: pd.DataFrame, hotspot_col: str, threshold_pct: float, nboot: int, rng: np.random.Generator) -> list[dict[str, object]]:
    hot_indices, ctrl_indices = [], []
    for _, q in d.groupby("stratum", sort=False):
        h = q.index[q[hotspot_col]].to_numpy(int)
        c = q.index[~q[hotspot_col]].to_numpy(int)
        if len(h) == 0 or len(c) == 0:
            continue
        for hi in h:
            selected = rng.choice(c, size=min(5, len(c)), replace=False)
            hot_indices.extend([hi] * len(selected))
            ctrl_indices.extend(selected.tolist())
    hot_indices = np.asarray(hot_indices, dtype=int)
    ctrl_indices = np.asarray(ctrl_indices, dtype=int)
    rows = []
    if len(hot_indices) == 0:
        return rows
    blocks = d.loc[hot_indices, "block"].astype(str).to_numpy()
    unique_blocks = np.array(sorted(np.unique(blocks)))
    for feature, label in PRIMARY_FEATURES.items():
        values = d[f"z_{feature}"].to_numpy(float)
        diff = values[hot_indices] - values[ctrl_indices]
        valid = np.isfinite(diff)
        diff, use_blocks = diff[valid], blocks[valid]
        block_sum = np.array([diff[use_blocks == b].sum() for b in unique_blocks])
        block_n = np.array([(use_blocks == b).sum() for b in unique_blocks])
        point = float(diff.mean())
        boots = np.empty(nboot)
        for i in range(nboot):
            sampled = rng.integers(0, len(unique_blocks), len(unique_blocks))
            boots[i] = block_sum[sampled].sum() / block_n[sampled].sum()
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append({
            "hotspot_threshold_percent": threshold_pct, "feature": feature, "feature_label": label,
            "n_unique_matched_hotspots": len(np.unique(hot_indices)), "n_matched_pairs": len(diff), "n_hotspot_blocks": len(unique_blocks),
            "standardized_hotspot_minus_control": point, "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(boots),
        })
    return rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    features = prepare_features(pd.read_csv(args.features, sep="\t"))
    matched_rows = []
    for pct, col in [(2.5, "hotspot_top_2_5pct"), (5.0, "hotspot_top_5pct"), (10.0, "hotspot_top_10pct")]:
        matched_rows.extend(matched_analysis(features, col, pct, args.bootstrap, rng))
    matched = pd.DataFrame(matched_rows)
    matched["fdr_within_threshold"] = matched.groupby("hotspot_threshold_percent")["bootstrap_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    matched.to_csv(outdir / "matched_hotspot_enrichment_012.tsv", sep="\t", index=False)

    link_rows = []
    ocr = pd.read_csv(args.ocr_metrics, sep="\t")
    same = np.sign(ocr["delta_rt_cpp8_1"]) == np.sign(ocr["delta_rt_cpp8_3"])
    ocr["consensus_instability"] = np.where(same, np.minimum(np.abs(ocr["delta_rt_cpp8_1"]), np.abs(ocr["delta_rt_cpp8_3"])), 0.0)
    ocr["hotspot"] = ocr["consensus_instability"] >= ocr["consensus_instability"].quantile(0.95)
    ocr["log_length"] = np.log1p(ocr["length"])
    for allele in ("cpp8_1", "cpp8_3"):
        for stage in ("ES", "MS", "LS"):
            sl = stage.lower()
            d = ocr[["chr", "block", "hotspot", "WT_wrt", "log_length", f"wt_{sl}_log_atac", f"delta_atac_{allele}_{sl}"]].copy()
            d.columns = ["chr", "block", "hotspot", "wt_wrt", "log_length", "wt_log_atac", "delta_atac"]
            d["abs_delta_atac"] = np.abs(d["delta_atac"])
            for col in ["abs_delta_atac", "wt_wrt", "log_length", "wt_log_atac"]:
                d[f"z_{col}"] = zscore(d[col])
            fit = smf.ols("z_abs_delta_atac ~ hotspot + z_wt_wrt + z_log_length + z_wt_log_atac + C(chr)", data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
            ci = fit.conf_int().loc["hotspot[T.True]"]
            link_rows.append({"layer": "OCR_abs_delta_ATAC", "allele": allele, "stage": stage, "promoter_window_bp": 0, "n_loci": len(d), "hotspot_standardized_difference": float(fit.params["hotspot[T.True]"]), "cluster_ci_low": float(ci.iloc[0]), "cluster_ci_high": float(ci.iloc[1]), "cluster_p": float(fit.pvalues["hotspot[T.True]"])})

    gm = pd.read_csv(args.gene_metrics, sep="\t")
    core = gm[(gm["stage"] == "MS") & (gm["promoter_window_bp"] == 1000)].copy()
    wide_rt = core.pivot_table(index=["gene_id", "chr", "tss", "block"], columns="allele", values="delta_rt").reset_index()
    same = np.sign(wide_rt["cpp8_1"]) == np.sign(wide_rt["cpp8_3"])
    wide_rt["consensus_instability"] = np.where(same, np.minimum(np.abs(wide_rt["cpp8_1"]), np.abs(wide_rt["cpp8_3"])), 0.0)
    wide_rt["hotspot"] = wide_rt["consensus_instability"] >= wide_rt["consensus_instability"].quantile(0.95)
    for allele in ("cpp8_1", "cpp8_3"):
        d = core[core["allele"] == allele].merge(wide_rt[["gene_id", "hotspot"]], on="gene_id")
        d["abs_rna_logfc"] = np.abs(d["rna_logfc"])
        for col in ["abs_rna_logfc", "baseline_rna", "wt_wrt", "wt_log_atac"]:
            d[f"z_{col}"] = zscore(d[col])
        fit = smf.ols("z_abs_rna_logfc ~ hotspot + z_baseline_rna + z_wt_wrt + z_wt_log_atac + C(chr)", data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
        ci = fit.conf_int().loc["hotspot[T.True]"]
        link_rows.append({"layer": "gene_abs_RNA_logFC", "allele": allele, "stage": "MS", "promoter_window_bp": 1000, "n_loci": len(d), "hotspot_standardized_difference": float(fit.params["hotspot[T.True]"]), "cluster_ci_low": float(ci.iloc[0]), "cluster_ci_high": float(ci.iloc[1]), "cluster_p": float(fit.pvalues["hotspot[T.True]"])})
    links = pd.DataFrame(link_rows)
    links["fdr_across_eight_links"] = multipletests(links["cluster_p"], method="fdr_bh")[1]
    links.to_csv(outdir / "hotspot_atac_rna_links_012.tsv", sep="\t", index=False)

    print("Matched hotspot enrichment")
    print(matched.to_string(index=False))
    print("\nATAC/RNA links")
    print(links.to_string(index=False))


if __name__ == "__main__":
    main()
