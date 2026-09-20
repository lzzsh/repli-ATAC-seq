#!/usr/bin/env python3
"""Fit the 3x3 WT RT-class by ATAC-stage coupling model with block bootstrap."""

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
RT_CLASSES = ("E", "M", "L")
STAGES = ("ES", "MS", "LS")
CELLS = tuple(f"{r}_{s}" for r in RT_CLASSES for s in STAGES)
FORMULA = (
    "z_delta_atac ~ 0 + C(cell) + C(cell):z_delta_rt + z_wt_wrt + "
    "z_wt_log_atac + z_rt_confidence + z_log_length + C(chr)"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--output-dir", required=True)
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


def slope_name(cell: str) -> str:
    return f"C(cell)[{cell}]:z_delta_rt"


def contrast_values(beta: dict[str, np.ndarray | float]) -> dict[str, np.ndarray | float]:
    diagonal = np.mean([beta["E_ES"], beta["M_MS"], beta["L_LS"]], axis=0)
    off = np.mean([beta[c] for c in CELLS if c not in {"E_ES", "M_MS", "L_LS"}], axis=0)
    lag = np.mean([beta["E_MS"], beta["M_LS"]], axis=0) - np.mean([beta["E_ES"], beta["M_MS"]], axis=0)
    ms_dom = np.mean([beta["E_MS"], beta["M_MS"], beta["L_MS"]], axis=0) - np.mean(
        [beta[f"{r}_{s}"] for r in RT_CLASSES for s in ("ES", "LS")], axis=0
    )
    return {"diagonal_dominance": diagonal - off, "post_replication_lag": lag, "global_ms_dominance": ms_dom}


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    data = pd.read_csv(args.matrix, sep="\t")
    data["cell"] = pd.Categorical(data["cell"], categories=CELLS)
    data["log_length"] = np.log1p(data["length"])
    prepared = {}
    for allele in ALLELES:
        d = data[data["allele"] == allele].copy()
        d["z_delta_rt"] = zscore(d["delta_rt"])
        d["z_wt_wrt"] = zscore(d["WT_wrt"])
        d["z_rt_confidence"] = zscore(d["rt_confidence"])
        d["z_log_length"] = zscore(d["log_length"])
        d["z_delta_atac"] = d.groupby("atac_stage", observed=True)["delta_atac"].transform(zscore)
        d["z_wt_log_atac"] = d.groupby("atac_stage", observed=True)["wt_log_atac"].transform(zscore)
        d = d.dropna(subset=["z_delta_atac", "z_delta_rt", "z_wt_wrt", "z_wt_log_atac", "z_rt_confidence", "z_log_length", "cell"]).copy()
        prepared[allele] = d

    blocks = np.array(sorted(set().union(*[set(d["block"].astype(str)) for d in prepared.values()])))
    fits = {}
    suffs = {}
    points = {}
    for allele in ALLELES:
        d = prepared[allele]
        fits[allele] = smf.ols(FORMULA, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
        suffs[allele] = sufficient(FORMULA, d, blocks)
        coef = solve(suffs[allele], np.arange(len(blocks)))
        points[allele] = {cell: float(coef[suffs[allele]["columns"].index(slope_name(cell))]) for cell in CELLS}

    boots = {allele: {cell: np.empty(args.bootstrap) for cell in CELLS} for allele in ALLELES}
    for i in range(args.bootstrap):
        sampled = rng.integers(0, len(blocks), len(blocks))
        for allele in ALLELES:
            coef = solve(suffs[allele], sampled)
            for cell in CELLS:
                boots[allele][cell][i] = coef[suffs[allele]["columns"].index(slope_name(cell))]

    cell_rows = []
    for allele in ALLELES:
        d = prepared[allele]
        for cell in CELLS:
            rt_class, stage = cell.split("_")
            values = boots[allele][cell]
            lo, hi = np.percentile(values, [2.5, 97.5])
            name = slope_name(cell)
            cell_rows.append({
                "allele": allele, "wt_rt_class": rt_class, "atac_stage": stage, "cell": cell,
                "n_ocr": int((d["cell"].astype(str) == cell).sum()), "n_blocks": len(blocks),
                "standardized_beta_delta_rt": points[allele][cell],
                "cluster_se": float(fits[allele].bse[name]), "cluster_p": float(fits[allele].pvalues[name]),
                "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(values),
            })
    cells = pd.DataFrame(cell_rows)
    cells["fdr_across_18_cells"] = multipletests(cells["cluster_p"], method="fdr_bh")[1]
    cells.to_csv(outdir / "stage_match_coefficients_010.tsv", sep="\t", index=False)

    contrast_rows = []
    contrast_boots = {}
    for allele in ALLELES:
        point_contrasts = contrast_values(points[allele])
        boot_contrasts = contrast_values(boots[allele])
        contrast_boots[allele] = boot_contrasts
        for name, value in point_contrasts.items():
            bv = np.asarray(boot_contrasts[name])
            lo, hi = np.percentile(bv, [2.5, 97.5])
            contrast_rows.append({"allele": allele, "contrast": name, "estimate": float(value), "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(bv)})
    contrasts = pd.DataFrame(contrast_rows)
    contrasts["fdr_across_six_allele_contrasts"] = multipletests(contrasts["bootstrap_p"], method="fdr_bh")[1]
    contrasts.to_csv(outdir / "stage_match_primary_contrasts_010.tsv", sep="\t", index=False)

    hetero_rows = []
    for name in contrast_boots[ALLELES[0]]:
        diff = np.asarray(contrast_boots["cpp8_1"][name]) - np.asarray(contrast_boots["cpp8_3"][name])
        lo, hi = np.percentile(diff, [2.5, 97.5])
        hetero_rows.append({
            "contrast": name,
            "cpp8_1_estimate": float(contrasts[(contrasts["allele"] == "cpp8_1") & (contrasts["contrast"] == name)]["estimate"].iloc[0]),
            "cpp8_3_estimate": float(contrasts[(contrasts["allele"] == "cpp8_3") & (contrasts["contrast"] == name)]["estimate"].iloc[0]),
            "difference_cpp8_1_minus_cpp8_3": float(np.mean(diff)), "difference_boot_ci_low": float(lo), "difference_boot_ci_high": float(hi), "difference_boot_p": bootstrap_p(diff),
        })
    hetero = pd.DataFrame(hetero_rows)
    hetero["fdr_across_three_contrasts"] = multipletests(hetero["difference_boot_p"], method="fdr_bh")[1]
    hetero.to_csv(outdir / "stage_match_allele_heterogeneity_010.tsv", sep="\t", index=False)

    lag_rows = []
    for allele in ALLELES:
        for lag in range(-2, 3):
            use = [f"{r}_{s}" for ri, r in enumerate(RT_CLASSES) for si, s in enumerate(STAGES) if si - ri == lag]
            point = float(np.mean([points[allele][cell] for cell in use]))
            bv = np.mean([boots[allele][cell] for cell in use], axis=0)
            lo, hi = np.percentile(bv, [2.5, 97.5])
            lag_rows.append({"allele": allele, "lag": lag, "n_cells": len(use), "cells": ",".join(use), "mean_beta": point, "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(bv)})
    pd.DataFrame(lag_rows).to_csv(outdir / "stage_lag_distance_010.tsv", sep="\t", index=False)

    print("Cell coefficients")
    print(cells.to_string(index=False))
    print("\nPrimary contrasts")
    print(contrasts.to_string(index=False))
    print("\nAllele heterogeneity")
    print(hetero.to_string(index=False))


if __name__ == "__main__":
    main()
