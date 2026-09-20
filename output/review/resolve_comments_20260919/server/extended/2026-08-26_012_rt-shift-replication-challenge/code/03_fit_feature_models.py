#!/usr/bin/env python3
"""Fit multivariable and signed-class models for consensus CPP8 RT instability."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import bs  # noqa: F401; exposed to formula evaluation
from statsmodels.stats.multitest import multipletests


PRIMARY_FEATURES = {
    "z_log_nearest_gene_length": "Nearest-gene length",
    "z_wt_expression_50kb_log2cpm": "WT expression within 50 kb",
    "z_class_i_te_fraction": "Class I TE fraction",
    "z_class_ii_te_fraction": "Class II TE fraction",
    "z_gc_fraction": "GC fraction",
    "z_log_distance_to_centromere_proxy": "Distance to centromere-repeat proxy",
    "z_log_gene_count_50kb": "Gene density within 50 kb",
    "z_log_ocr_count_50kb": "OCR density within 50 kb",
}
SECONDARY_FEATURES = {
    "z_h3k27ac_fraction": "H3K27ac fraction",
    "z_h3k4me1_fraction": "H3K4me1 fraction",
    "z_h3k4me3_fraction": "H3K4me3 fraction",
    "z_h3k9me2_fraction": "H3K9me2 fraction",
    "z_h3k27me3_fraction": "H3K27me3 fraction",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def prepare(data: pd.DataFrame) -> pd.DataFrame:
    d = data.copy()
    d["log_nearest_gene_length"] = np.log1p(d["nearest_gene_length"])
    d["log_distance_to_centromere_proxy"] = np.log1p(d["distance_to_centromere_proxy"])
    d["log_gene_count_50kb"] = np.log1p(d["gene_count_50kb"])
    d["log_ocr_count_50kb"] = np.log1p(d["ocr_count_50kb"])
    d["z_consensus_instability"] = zscore(d["consensus_instability"])
    d["z_log_wt_g1_coverage"] = zscore(d["log_wt_g1_coverage"])
    source_cols = {
        "z_log_nearest_gene_length": "log_nearest_gene_length",
        "z_wt_expression_50kb_log2cpm": "wt_expression_50kb_log2cpm",
        "z_class_i_te_fraction": "class_i_te_fraction",
        "z_class_ii_te_fraction": "class_ii_te_fraction",
        "z_gc_fraction": "gc_fraction",
        "z_log_distance_to_centromere_proxy": "log_distance_to_centromere_proxy",
        "z_log_gene_count_50kb": "log_gene_count_50kb",
        "z_log_ocr_count_50kb": "log_ocr_count_50kb",
        "z_h3k27ac_fraction": "h3k27ac_fraction",
        "z_h3k4me1_fraction": "h3k4me1_fraction",
        "z_h3k4me3_fraction": "h3k4me3_fraction",
        "z_h3k9me2_fraction": "h3k9me2_fraction",
        "z_h3k27me3_fraction": "h3k27me3_fraction",
    }
    for zcol, source in source_cols.items():
        d[zcol] = zscore(d[source])
    return d


def coefficient_rows(fit, features: dict[str, str], family: str, n: int, n_blocks: int) -> list[dict[str, object]]:
    rows = []
    ci = fit.conf_int()
    for col, label in features.items():
        rows.append({
            "family": family, "feature": col.removeprefix("z_"), "feature_label": label,
            "n_bins": n, "n_blocks": n_blocks, "standardized_beta": float(fit.params[col]),
            "cluster_se": float(fit.bse[col]), "cluster_ci_low": float(ci.loc[col, 0]), "cluster_ci_high": float(ci.loc[col, 1]),
            "cluster_p": float(fit.pvalues[col]), "adjusted_r2": float(fit.rsquared_adj),
        })
    return rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    d = prepare(pd.read_csv(args.features, sep="\t"))

    core = " + ".join(PRIMARY_FEATURES)
    formula = f"z_consensus_instability ~ bs(WT_wrt, df=4) + z_log_wt_g1_coverage + {core} + C(chr)"
    use = d.dropna(subset=["z_consensus_instability", "WT_wrt", "z_log_wt_g1_coverage", *PRIMARY_FEATURES]).copy()
    fit = smf.ols(formula, data=use).fit(cov_type="cluster", cov_kwds={"groups": use["block"]})
    rows = coefficient_rows(fit, PRIMARY_FEATURES, "primary", len(use), use["block"].nunique())

    secondary = " + ".join(SECONDARY_FEATURES)
    sec_formula = f"z_consensus_instability ~ bs(WT_wrt, df=4) + z_log_wt_g1_coverage + {core} + {secondary} + C(chr)"
    sec_use = d.dropna(subset=["z_consensus_instability", "WT_wrt", "z_log_wt_g1_coverage", *PRIMARY_FEATURES, *SECONDARY_FEATURES]).copy()
    sec_fit = smf.ols(sec_formula, data=sec_use).fit(cov_type="cluster", cov_kwds={"groups": sec_use["block"]})
    rows.extend(coefficient_rows(sec_fit, SECONDARY_FEATURES, "secondary_histone", len(sec_use), sec_use["block"].nunique()))
    results = pd.DataFrame(rows)
    results["fdr_within_family"] = results.groupby("family")["cluster_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    results.to_csv(outdir / "feature_multivariable_models_012.tsv", sep="\t", index=False)

    class_order = ["stable", "concordant_earlier", "concordant_later", "discordant"]
    d["shift_class"] = pd.Categorical(d["signed_shift_class"], categories=class_order)
    signed_rows = []
    all_features = {**PRIMARY_FEATURES, **SECONDARY_FEATURES}
    for feature, label in all_features.items():
        q = d.dropna(subset=[feature, "shift_class", "WT_wrt", "z_log_wt_g1_coverage"]).copy()
        f = smf.ols(f"{feature} ~ C(shift_class) + bs(WT_wrt, df=4) + z_log_wt_g1_coverage + C(chr)", data=q).fit(cov_type="cluster", cov_kwds={"groups": q["block"]})
        ci = f.conf_int()
        for cls in class_order[1:]:
            name = f"C(shift_class)[T.{cls}]"
            signed_rows.append({
                "feature": feature.removeprefix("z_"), "feature_label": label, "shift_class_vs_stable": cls,
                "n_bins": len(q), "standardized_difference": float(f.params[name]), "cluster_se": float(f.bse[name]),
                "cluster_ci_low": float(ci.loc[name, 0]), "cluster_ci_high": float(ci.loc[name, 1]), "cluster_p": float(f.pvalues[name]),
            })
    signed = pd.DataFrame(signed_rows)
    signed["fdr_across_feature_by_class_tests"] = multipletests(signed["cluster_p"], method="fdr_bh")[1]
    signed.to_csv(outdir / "signed_shift_class_models_012.tsv", sep="\t", index=False)

    print("Multivariable feature model")
    print(results.to_string(index=False))
    print("\nSigned shift classes")
    print(signed.to_string(index=False))


if __name__ == "__main__":
    main()
