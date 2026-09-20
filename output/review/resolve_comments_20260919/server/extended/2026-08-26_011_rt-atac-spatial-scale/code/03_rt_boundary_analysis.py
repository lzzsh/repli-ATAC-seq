#!/usr/bin/env python3
"""Test whether local RT-ATAC coupling differs near candidate WT RT boundaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrices
from statsmodels.stats.multitest import multipletests


SEED = 20260826
ALLELES = ("cpp8_1", "cpp8_3")
STAGES = ("ES", "MS", "LS")
TOP_PCT = (5, 10, 15)
FORMULA = "z_delta_atac ~ z_delta_rt * C(boundary_group) + z_wt_wrt + z_wt_log_atac + z_log_length + C(chr)"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--metrics", required=True)
    p.add_argument("--rt-bins", required=True)
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


def sufficient(data: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(FORMULA, data, return_type="dataframe", NA_action="drop")
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


def nearest_distance(query: np.ndarray, positions: np.ndarray) -> np.ndarray:
    if len(positions) == 0:
        return np.full(len(query), np.nan)
    idx = np.searchsorted(positions, query)
    left = positions[np.maximum(idx - 1, 0)]
    right = positions[np.minimum(idx, len(positions) - 1)]
    return np.minimum(np.abs(query - left), np.abs(query - right))


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.metrics, sep="\t")
    rt = pd.read_csv(args.rt_bins, sep="\t")

    boundary_candidates = []
    for chrom, d in rt.groupby("chr", sort=False):
        d = d.sort_values("start").copy()
        d["wt_wrt_20kb_rolling_median"] = d["WT_wrt"].rolling(20, center=True, min_periods=10).median()
        d["abs_first_difference"] = d["wt_wrt_20kb_rolling_median"].diff().abs()
        # A threshold alone would label long runs of adjacent 1-kb bins as many
        # separate boundaries.  Collapse each smoothed transition to its local
        # maximum within a 20-kb neighbourhood before applying the percentile.
        local_max = d["abs_first_difference"].rolling(21, center=True, min_periods=1).max()
        d["is_local_max_20kb"] = d["abs_first_difference"].eq(local_max)
        for row in d.dropna(subset=["abs_first_difference"])[["start", "abs_first_difference", "is_local_max_20kb"]].itertuples(index=False):
            boundary_candidates.append({"chr": chrom, "position": int(row.start), "abs_first_difference": float(row.abs_first_difference), "is_local_max_20kb": bool(row.is_local_max_20kb)})
    candidates = pd.DataFrame(boundary_candidates)
    candidates.to_csv(outdir / "candidate_rt_boundaries_011.tsv.gz", sep="\t", index=False)

    rows = []
    blocks = np.array(sorted(metrics["block"].astype(str).unique()))
    for top_pct in TOP_PCT:
        threshold = float(np.percentile(candidates["abs_first_difference"], 100 - top_pct))
        selected = candidates[(candidates["abs_first_difference"] >= threshold) & candidates["is_local_max_20kb"]]
        distance = np.full(len(metrics), np.nan)
        for chrom, oi in metrics.groupby("chr", sort=False).groups.items():
            positions = np.sort(selected.loc[selected["chr"] == chrom, "position"].to_numpy(int))
            distance[np.asarray(oi, dtype=int)] = nearest_distance(metrics.loc[oi, "midpoint"].to_numpy(int), positions)
        boundary_group = np.where(distance <= 10_000, "near", np.where(distance >= 50_000, "interior", "excluded"))
        for allele in ALLELES:
            for stage in STAGES:
                sl = stage.lower()
                d = metrics[["chr", "block", "length", "WT_wrt", f"delta_rt_{allele}", f"wt_{sl}_log_atac", f"delta_atac_{allele}_{sl}"]].copy()
                d.columns = ["chr", "block", "length", "wt_wrt", "delta_rt", "wt_log_atac", "delta_atac"]
                d["boundary_group"] = pd.Categorical(boundary_group, categories=["interior", "near", "excluded"])
                d = d[d["boundary_group"].isin(["interior", "near"])].copy()
                d["boundary_group"] = d["boundary_group"].cat.remove_unused_categories()
                d["log_length"] = np.log1p(d["length"])
                d = d.dropna().copy()
                for col in ["delta_atac", "delta_rt", "wt_wrt", "wt_log_atac", "log_length"]:
                    d[f"z_{col}"] = zscore(d[col])
                fit = smf.ols(FORMULA, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
                suff = sufficient(d, blocks)
                base_name = "z_delta_rt"
                int_name = "z_delta_rt:C(boundary_group)[T.near]"
                point_coef = solve(suff, np.arange(len(blocks)))
                base_idx = suff["columns"].index(base_name)
                int_idx = suff["columns"].index(int_name)
                point_interior = float(point_coef[base_idx])
                point_diff = float(point_coef[int_idx])
                boot_interior, boot_diff = np.empty(args.bootstrap), np.empty(args.bootstrap)
                for i in range(args.bootstrap):
                    b = solve(suff, rng.integers(0, len(blocks), len(blocks)))
                    boot_interior[i], boot_diff[i] = b[base_idx], b[int_idx]
                boot_near = boot_interior + boot_diff
                for group, point, bv in [("interior", point_interior, boot_interior), ("near", point_interior + point_diff, boot_near), ("near_minus_interior", point_diff, boot_diff)]:
                    lo, hi = np.percentile(bv, [2.5, 97.5])
                    rows.append({
                        "boundary_top_percent": top_pct, "boundary_threshold": threshold, "n_boundaries": len(selected),
                        "allele": allele, "stage": stage, "effect": group,
                        "n_ocr_near": int((d["boundary_group"] == "near").sum()), "n_ocr_interior": int((d["boundary_group"] == "interior").sum()),
                        "estimate": point, "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(bv),
                        "cluster_interaction_p": float(fit.pvalues[int_name]) if group == "near_minus_interior" else np.nan,
                    })
    result = pd.DataFrame(rows)
    interaction = result["effect"] == "near_minus_interior"
    result.loc[interaction, "interaction_fdr_across_18_tests"] = multipletests(result.loc[interaction, "cluster_interaction_p"], method="fdr_bh")[1]
    result.to_csv(outdir / "rt_boundary_sensitivity_011.tsv", sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
