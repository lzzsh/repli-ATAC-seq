#!/usr/bin/env python3
"""Fit nested-scale, annulus and offset CPP8 RT-ATAC association models."""

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
SCALES = (1, 5, 25, 100)
OFFSETS = tuple(range(-100, 101, 5))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--metrics", required=True)
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


def trend_weights(x: np.ndarray) -> np.ndarray:
    centered = x - x.mean()
    return centered / np.sum(centered**2)


def offset_tag(offset: int) -> str:
    return f"m{abs(offset)}" if offset < 0 else f"p{offset}"


def prepare_scale_data(metrics: pd.DataFrame, allele: str, stage: str, scale: int) -> pd.DataFrame:
    sl = stage.lower()
    d = metrics[["chr", "block", "length", f"delta_atac_{allele}_{sl}", f"wt_{sl}_log_atac", f"delta_rt_{allele}_scale_{scale}kb", f"WT_wrt_scale_{scale}kb"]].copy()
    d.columns = ["chr", "block", "length", "delta_atac", "wt_log_atac", "delta_rt", "wt_wrt"]
    d["log_length"] = np.log1p(d["length"])
    d = d.dropna().copy()
    for col in ["delta_atac", "wt_log_atac", "delta_rt", "wt_wrt", "log_length"]:
        d[f"z_{col}"] = zscore(d[col])
    return d


SCALE_FORMULA = "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + z_log_length + C(chr)"


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.metrics, sep="\t")
    blocks = np.array(sorted(metrics["block"].astype(str).unique()))

    coefficient_rows = []
    contrast_rows = []
    all_boots: dict[tuple[str, str], dict[int, np.ndarray]] = {}
    for allele in ALLELES:
        for stage in STAGES:
            suffs, points, fits = {}, {}, {}
            for scale in SCALES:
                d = prepare_scale_data(metrics, allele, stage, scale)
                fits[scale] = smf.ols(SCALE_FORMULA, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
                suffs[scale] = sufficient(SCALE_FORMULA, d, blocks)
                coef = solve(suffs[scale], np.arange(len(blocks)))
                points[scale] = float(coef[suffs[scale]["columns"].index("z_delta_rt")])
            boots = {scale: np.empty(args.bootstrap) for scale in SCALES}
            for i in range(args.bootstrap):
                sampled = rng.integers(0, len(blocks), len(blocks))
                for scale in SCALES:
                    coef = solve(suffs[scale], sampled)
                    boots[scale][i] = coef[suffs[scale]["columns"].index("z_delta_rt")]
            all_boots[(allele, stage)] = boots
            for scale in SCALES:
                lo, hi = np.percentile(boots[scale], [2.5, 97.5])
                coefficient_rows.append({
                    "allele": allele, "stage": stage, "scale_radius_kb": scale,
                    "n_ocr": int(fits[scale].nobs), "n_blocks": len(blocks),
                    "standardized_beta_delta_rt": points[scale], "cluster_se": float(fits[scale].bse["z_delta_rt"]),
                    "cluster_p": float(fits[scale].pvalues["z_delta_rt"]), "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(boots[scale]),
                    "adjusted_r2": float(fits[scale].rsquared_adj),
                })
            local = points[1] - np.mean([points[25], points[100]])
            local_boot = boots[1] - (boots[25] + boots[100]) / 2
            llo, lhi = np.percentile(local_boot, [2.5, 97.5])
            x = np.log10(np.asarray(SCALES, dtype=float))
            weights = trend_weights(x)
            trend = float(np.dot(weights, np.asarray([points[s] for s in SCALES])))
            trend_boot = sum(weights[i] * boots[s] for i, s in enumerate(SCALES))
            tlo, thi = np.percentile(trend_boot, [2.5, 97.5])
            contrast_rows.extend([
                {"allele": allele, "stage": stage, "contrast": "local_dominance", "estimate": local, "bootstrap_ci_low": float(llo), "bootstrap_ci_high": float(lhi), "bootstrap_p": bootstrap_p(local_boot)},
                {"allele": allele, "stage": stage, "contrast": "beta_trend_per_log10_radius", "estimate": trend, "bootstrap_ci_low": float(tlo), "bootstrap_ci_high": float(thi), "bootstrap_p": bootstrap_p(trend_boot)},
            ])

    coefficients = pd.DataFrame(coefficient_rows)
    coefficients["fdr_across_24_scale_tests"] = multipletests(coefficients["cluster_p"], method="fdr_bh")[1]
    coefficients.to_csv(outdir / "multiscale_rt_atac_coefficients_011.tsv", sep="\t", index=False)
    contrasts = pd.DataFrame(contrast_rows)
    contrasts["fdr_within_contrast_family"] = contrasts.groupby("contrast")["bootstrap_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    contrasts.to_csv(outdir / "multiscale_primary_contrasts_011.tsv", sep="\t", index=False)

    annulus_rows = []
    annulus_names = ["local", "annulus_1_5kb", "annulus_5_20kb", "annulus_20_50kb", "annulus_50_100kb"]
    annulus_formula = "z_delta_atac ~ z_local + z_annulus_1_5kb + z_annulus_5_20kb + z_annulus_20_50kb + z_annulus_50_100kb + z_wt_wrt + z_wt_log_atac + z_log_length + C(chr)"
    for allele in ALLELES:
        for stage in STAGES:
            sl = stage.lower()
            d = metrics[["chr", "block", "length", f"delta_atac_{allele}_{sl}", f"wt_{sl}_log_atac", f"delta_rt_{allele}", f"delta_rt_{allele}_annulus_1_5kb", f"delta_rt_{allele}_annulus_5_20kb", f"delta_rt_{allele}_annulus_20_50kb", f"delta_rt_{allele}_annulus_50_100kb", "WT_wrt_scale_100kb"]].copy()
            d.columns = ["chr", "block", "length", "delta_atac", "wt_log_atac", "local", "annulus_1_5kb", "annulus_5_20kb", "annulus_20_50kb", "annulus_50_100kb", "wt_wrt"]
            d["log_length"] = np.log1p(d["length"])
            d = d.dropna().copy()
            for col in ["delta_atac", "wt_log_atac", "local", "annulus_1_5kb", "annulus_5_20kb", "annulus_20_50kb", "annulus_50_100kb", "wt_wrt", "log_length"]:
                d[f"z_{col}"] = zscore(d[col])
            fit = smf.ols(annulus_formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
            suff = sufficient(annulus_formula, d, blocks)
            coef = solve(suff, np.arange(len(blocks)))
            boot = {name: np.empty(args.bootstrap) for name in annulus_names}
            for i in range(args.bootstrap):
                b = solve(suff, rng.integers(0, len(blocks), len(blocks)))
                for name in annulus_names:
                    boot[name][i] = b[suff["columns"].index(f"z_{name}")]
            for name in annulus_names:
                point = float(coef[suff["columns"].index(f"z_{name}")])
                lo, hi = np.percentile(boot[name], [2.5, 97.5])
                annulus_rows.append({"allele": allele, "stage": stage, "spatial_component": name, "n_ocr": len(d), "standardized_beta": point, "cluster_se": float(fit.bse[f"z_{name}"]), "cluster_p": float(fit.pvalues[f"z_{name}"]), "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi), "bootstrap_p": bootstrap_p(boot[name])})
    annulus = pd.DataFrame(annulus_rows)
    annulus["fdr_across_30_components"] = multipletests(annulus["cluster_p"], method="fdr_bh")[1]
    annulus.to_csv(outdir / "local_vs_annulus_models_011.tsv", sep="\t", index=False)

    offset_rows = []
    for allele in ALLELES:
        for stage in STAGES:
            sl = stage.lower()
            for offset in OFFSETS:
                tag = offset_tag(offset)
                d = metrics[["chr", "block", "length", f"delta_atac_{allele}_{sl}", f"wt_{sl}_log_atac", f"delta_rt_{allele}_offset_{tag}kb", f"WT_wrt_offset_{tag}kb"]].copy()
                d.columns = ["chr", "block", "length", "delta_atac", "wt_log_atac", "delta_rt", "wt_wrt"]
                d["log_length"] = np.log1p(d["length"])
                d = d.dropna().copy()
                for col in ["delta_atac", "wt_log_atac", "delta_rt", "wt_wrt", "log_length"]:
                    d[f"z_{col}"] = zscore(d[col])
                fit = smf.ols(SCALE_FORMULA, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
                ci = fit.conf_int().loc["z_delta_rt"]
                offset_rows.append({"allele": allele, "stage": stage, "offset_kb": offset, "window_width_kb": 5, "n_ocr": int(fit.nobs), "standardized_beta_delta_rt": float(fit.params["z_delta_rt"]), "cluster_ci_low": float(ci.iloc[0]), "cluster_ci_high": float(ci.iloc[1]), "cluster_p": float(fit.pvalues["z_delta_rt"])})
    pd.DataFrame(offset_rows).to_csv(outdir / "spatial_offset_curve_011.tsv.gz", sep="\t", index=False)

    print("Nested scale coefficients")
    print(coefficients.to_string(index=False))
    print("\nPrimary contrasts")
    print(contrasts.to_string(index=False))
    print("\nSimultaneous annulus components")
    print(annulus.to_string(index=False))


if __name__ == "__main__":
    main()
