#!/usr/bin/env python3
"""Test spatially robust associations between mutant RT shifts and phased ATAC shifts.

The unit of analysis is an ATAC OCR. OCRs are overlap-weighted to 1-kb Repli-seq
bins. ATAC count columns are normalized with a DESeq-style positive-count median
ratio size factor. All inference clusters observations by 1-Mb genomic block.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrices
from statsmodels.stats.multitest import multipletests


SEED = 20260826
STAGES = ("ES", "MS", "LS")
MUTANTS = ("cpp8", "cpp11")
COLORS = {"cpp8": "#C25759", "cpp11": "#397CA3"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--atac", required=True)
    p.add_argument("--repli", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=1000)
    return p.parse_args()


def wrt(es: pd.Series, ms: pd.Series, ls: pd.Series) -> np.ndarray:
    esv = pd.to_numeric(es, errors="coerce").to_numpy(float)
    msv = pd.to_numeric(ms, errors="coerce").to_numpy(float)
    lsv = pd.to_numeric(ls, errors="coerce").to_numpy(float)
    denom = esv + msv + lsv
    return np.divide(0.5 * msv + lsv, denom, out=np.full_like(denom, np.nan), where=denom > 0)


def positive_count_size_factors(counts: np.ndarray) -> np.ndarray:
    """DESeq2-like median ratios using the positive-count geometric mean."""
    positive = counts > 0
    npos = positive.sum(axis=1)
    log_counts = np.zeros_like(counts, dtype=float)
    np.log(counts, out=log_counts, where=positive)
    geo = np.exp(log_counts.sum(axis=1) / np.maximum(npos, 1))
    usable = (npos >= 2) & np.isfinite(geo) & (geo > 0)
    sf = np.empty(counts.shape[1], dtype=float)
    for j in range(counts.shape[1]):
        valid = usable & positive[:, j]
        sf[j] = np.median(counts[valid, j] / geo[valid])
    sf = sf / np.exp(np.mean(np.log(sf)))
    return sf


def map_ocr_to_rt(ocr: pd.DataFrame, rt: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Overlap-weight RT-bin values for each OCR without requiring interval packages."""
    result = np.full((len(ocr), len(value_cols)), np.nan, dtype=float)
    overlap_bp = np.zeros(len(ocr), dtype=int)
    n_bins = np.zeros(len(ocr), dtype=int)

    for chrom, oi in ocr.groupby("chr", sort=False).groups.items():
        od = ocr.loc[oi].sort_values("start")
        rd = rt[rt["chr"] == chrom].sort_values("start")
        if rd.empty:
            continue
        rs = rd["start"].to_numpy(int)
        re = rd["end"].to_numpy(int)
        rv = rd[value_cols].to_numpy(float)
        for idx, start, end in od[["start", "end"]].itertuples():
            left = int(np.searchsorted(re, start, side="right"))
            right = int(np.searchsorted(rs, end, side="left"))
            if right <= left:
                continue
            ov = np.minimum(re[left:right], end) - np.maximum(rs[left:right], start)
            keep = ov > 0
            if not np.any(keep):
                continue
            weights = ov[keep].astype(float)
            vals = rv[left:right][keep]
            valid = np.isfinite(vals)
            denom = (valid * weights[:, None]).sum(axis=0)
            num = np.where(valid, vals, 0.0) * weights[:, None]
            mapped = np.divide(num.sum(axis=0), denom, out=np.full(len(value_cols), np.nan), where=denom > 0)
            result[idx, :] = mapped
            overlap_bp[idx] = int(weights.sum())
            n_bins[idx] = int(keep.sum())

    mapped = pd.DataFrame(result, columns=value_cols, index=ocr.index)
    mapped["rt_overlap_bp"] = overlap_bp
    mapped["n_rt_bins"] = n_bins
    return mapped


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def block_bootstrap_coefficient(
    formula: str,
    data: pd.DataFrame,
    coefficient: str,
    nboot: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    y, x = dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    groups = data.loc[x.index, "block"].astype(str)
    unique_groups, codes = np.unique(groups, return_inverse=True)
    p = x.shape[1]
    xtx = np.zeros((len(unique_groups), p, p), dtype=float)
    xty = np.zeros((len(unique_groups), p), dtype=float)
    xa = x.to_numpy(float)
    ya = y.iloc[:, 0].to_numpy(float)
    for g in range(len(unique_groups)):
        take = codes == g
        xg = xa[take]
        yg = ya[take]
        xtx[g] = xg.T @ xg
        xty[g] = xg.T @ yg
    coef_idx = x.columns.get_loc(coefficient)
    point_vec = np.linalg.pinv(xtx.sum(axis=0)) @ xty.sum(axis=0)
    point = float(point_vec[coef_idx])
    boots = np.empty(nboot, dtype=float)
    ng = len(unique_groups)
    for i in range(nboot):
        sampled = rng.integers(0, ng, ng)
        beta = np.linalg.pinv(xtx[sampled].sum(axis=0)) @ xty[sampled].sum(axis=0)
        boots[i] = beta[coef_idx]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def clustered_bin_summary(d: pd.DataFrame) -> pd.DataFrame:
    out = []
    use = d.copy()
    use["rt_decile"] = pd.qcut(use["delta_rt"], 10, labels=False, duplicates="drop") + 1
    for q, qd in use.groupby("rt_decile"):
        fit = smf.ols("delta_atac ~ 1", data=qd).fit(cov_type="cluster", cov_kwds={"groups": qd["block"]})
        out.append(
            {
                "rt_decile": int(q),
                "mean_delta_rt": float(qd["delta_rt"].mean()),
                "mean_delta_atac": float(fit.params["Intercept"]),
                "ci_low": float(fit.conf_int().loc["Intercept", 0]),
                "ci_high": float(fit.conf_int().loc["Intercept", 1]),
                "n_ocr": len(qd),
            }
        )
    return pd.DataFrame(out)


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    atac_names = [
        "chr", "start", "end",
        "wt_es_1", "wt_es_2", "wt_ms", "wt_ls",
        "cpp11_es_1", "cpp11_es_2", "cpp11_ms", "cpp11_ls",
        "cpp8_es_1", "cpp8_es_2", "cpp8_ms", "cpp8_ls",
    ]
    atac = pd.read_csv(args.atac, sep="\t", header=None, names=atac_names)
    atac = atac[atac["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy().reset_index(drop=True)
    atac[["start", "end"]] = atac[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    atac = atac.dropna(subset=["start", "end"]).reset_index(drop=True)
    atac[["start", "end"]] = atac[["start", "end"]].astype(int)
    sample_cols = atac_names[3:]
    counts = atac[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
    size_factors = positive_count_size_factors(counts)
    normalized = counts / size_factors[None, :]
    atac[sample_cols] = normalized
    pd.DataFrame({"sample": sample_cols, "size_factor": size_factors}).to_csv(
        outdir / "atac_positive_count_size_factors_004.tsv", sep="\t", index=False
    )

    repli_names = [
        "chr", "start", "end",
        "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
        "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
        "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
        "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
        "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS",
    ]
    raw_rt = pd.read_csv(args.repli, sep="\t", header=None, names=repli_names)
    raw_rt = raw_rt[raw_rt["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy()
    raw_rt[["start", "end"]] = raw_rt[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    raw_rt = raw_rt.dropna(subset=["start", "end"])
    raw_rt[["start", "end"]] = raw_rt[["start", "end"]].astype(int)
    library_cols = repli_names[3:]
    raw_rt[library_cols] = raw_rt[library_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    length_kb = (raw_rt["end"] - raw_rt["start"] + 1).to_numpy(float) / 1000.0
    rpk = raw_rt[library_cols].to_numpy(float) / length_kb[:, None]
    tpm = rpk / (rpk.sum(axis=0, keepdims=True) / 1_000_000.0)
    raw_rt[library_cols] = tpm

    rt = raw_rt[["chr", "start", "end"]].copy()
    replicate_map = {
        "WT": {"G1": ["WT-1-G1", "WT-2-G1"], "ES": ["WT-1-ES", "WT-2-ES"], "MS": ["WT-1-MS"], "LS": ["WT-1-LS"]},
        "sol1_8": {"G1": ["sol1_8-1-G1", "sol1_8-2-G1"], "ES": ["sol1_8-1-ES", "sol1_8-2-ES"], "MS": ["sol1_8-1-MS"], "LS": ["sol1_8-1-LS"]},
        "tcx2_3": {"G1": ["tcx2_3-1-G1"], "ES": ["tcx2_3-1-ES"], "MS": ["tcx2_3-1-MS"], "LS": ["tcx2_3-1-LS"]},
    }
    for sample, stage_map in replicate_map.items():
        g1 = raw_rt[stage_map["G1"]].mean(axis=1)
        for stage in STAGES:
            signal = raw_rt[stage_map[stage]].mean(axis=1)
            rt[f"{sample}_{stage}_norm"] = signal / (g1 + 1e-6)
    rt["wt_wrt"] = wrt(rt.WT_ES_norm, rt.WT_MS_norm, rt.WT_LS_norm)
    rt["cpp8_wrt"] = wrt(rt.tcx2_3_ES_norm, rt.tcx2_3_MS_norm, rt.tcx2_3_LS_norm)
    rt["cpp11_wrt"] = wrt(rt.sol1_8_ES_norm, rt.sol1_8_MS_norm, rt.sol1_8_LS_norm)
    rt["delta_rt_cpp8"] = rt["cpp8_wrt"] - rt["wt_wrt"]
    rt["delta_rt_cpp11"] = rt["cpp11_wrt"] - rt["wt_wrt"]

    mapped = map_ocr_to_rt(
        atac[["chr", "start", "end"]],
        rt[["chr", "start", "end", "wt_wrt", "delta_rt_cpp8", "delta_rt_cpp11"]],
        ["wt_wrt", "delta_rt_cpp8", "delta_rt_cpp11"],
    )
    metrics = pd.concat([atac[["chr", "start", "end"]], mapped], axis=1)
    metrics["length"] = metrics["end"] - metrics["start"]
    metrics["block"] = metrics["chr"].astype(str) + ":" + (metrics["start"] // 1_000_000).astype(str)

    group_cols = {
        "ES": {"wt": ["wt_es_1", "wt_es_2"], "cpp8": ["cpp8_es_1", "cpp8_es_2"], "cpp11": ["cpp11_es_1", "cpp11_es_2"]},
        "MS": {"wt": ["wt_ms"], "cpp8": ["cpp8_ms"], "cpp11": ["cpp11_ms"]},
        "LS": {"wt": ["wt_ls"], "cpp8": ["cpp8_ls"], "cpp11": ["cpp11_ls"]},
    }
    pseudocount = 0.5
    for stage, groups in group_cols.items():
        stage_l = stage.lower()
        for genotype, cols in groups.items():
            metrics[f"{genotype}_{stage_l}_atac"] = atac[cols].mean(axis=1)
        metrics[f"wt_{stage_l}_log_atac"] = np.log2(metrics[f"wt_{stage_l}_atac"] + pseudocount)
        for mutant in MUTANTS:
            metrics[f"delta_atac_{mutant}_{stage_l}"] = np.log2(
                metrics[f"{mutant}_{stage_l}_atac"] + pseudocount
            ) - np.log2(metrics[f"wt_{stage_l}_atac"] + pseudocount)

    valid_map = metrics["wt_wrt"].notna() & (metrics["rt_overlap_bp"] > 0)
    metrics = metrics.loc[valid_map].reset_index(drop=True)
    metrics.to_csv(outdir / "ocr_rt_atac_metrics_004.tsv.gz", sep="\t", index=False)

    formula = (
        "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + "
        "z_log_length + C(chr)"
    )
    rows = []
    binned = []
    stacked = []
    for mutant in MUTANTS:
        for stage in STAGES:
            stage_l = stage.lower()
            d = metrics[
                ["chr", "block", "length", "wt_wrt", f"delta_rt_{mutant}",
                 f"wt_{stage_l}_log_atac", f"delta_atac_{mutant}_{stage_l}"]
            ].copy()
            d = d.rename(
                columns={
                    f"delta_rt_{mutant}": "delta_rt",
                    f"wt_{stage_l}_log_atac": "wt_log_atac",
                    f"delta_atac_{mutant}_{stage_l}": "delta_atac",
                }
            ).dropna()
            d["log_length"] = np.log1p(d["length"])
            for col in ["delta_atac", "delta_rt", "wt_wrt", "wt_log_atac", "log_length"]:
                d[f"z_{col}"] = zscore(d[col])
            fit = smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
            boot_point, boot_lo, boot_hi = block_bootstrap_coefficient(
                formula, d, "z_delta_rt", args.bootstrap, rng
            )
            ci = fit.conf_int().loc["z_delta_rt"]
            rows.append(
                {
                    "mutant": mutant,
                    "stage": stage,
                    "n_ocr": len(d),
                    "n_blocks": d["block"].nunique(),
                    "standardized_beta_delta_rt": float(fit.params["z_delta_rt"]),
                    "cluster_se": float(fit.bse["z_delta_rt"]),
                    "cluster_ci_low": float(ci.iloc[0]),
                    "cluster_ci_high": float(ci.iloc[1]),
                    "cluster_p": float(fit.pvalues["z_delta_rt"]),
                    "block_bootstrap_beta": boot_point,
                    "bootstrap_ci_low": boot_lo,
                    "bootstrap_ci_high": boot_hi,
                    "adjusted_r2": float(fit.rsquared_adj),
                }
            )
            bs = clustered_bin_summary(d[["delta_rt", "delta_atac", "block"]])
            bs["mutant"] = mutant
            bs["stage"] = stage
            binned.append(bs)
            ds = d.copy()
            ds["mutant"] = mutant
            ds["stage"] = stage
            ds["stage_index"] = STAGES.index(stage)
            stacked.append(ds)

    assoc = pd.DataFrame(rows)
    assoc["fdr_across_six_tests"] = multipletests(assoc["cluster_p"], method="fdr_bh")[1]
    assoc.to_csv(outdir / "rt_atac_stage_associations_004.tsv", sep="\t", index=False)
    binned_df = pd.concat(binned, ignore_index=True)
    binned_df.to_csv(outdir / "rt_atac_decile_summaries_004.tsv", sep="\t", index=False)

    trend_rows = []
    # C(stage) supplies the stage main effects. Adding stage_index itself would
    # be exactly collinear with those dummies, so only its interaction with RT
    # is included to test a linear ES-to-LS change in the RT slope.
    trend_formula = (
        "z_delta_atac ~ z_delta_rt + z_delta_rt:stage_index + z_wt_wrt + "
        "z_wt_log_atac + z_log_length + C(chr) + C(stage)"
    )
    for mutant in MUTANTS:
        ds = pd.concat([x for x in stacked if x["mutant"].iloc[0] == mutant], ignore_index=True)
        fit = smf.ols(trend_formula, data=ds).fit(cov_type="cluster", cov_kwds={"groups": ds["block"]})
        coef_name = "z_delta_rt:stage_index"
        ci = fit.conf_int().loc[coef_name]
        point, lo, hi = block_bootstrap_coefficient(trend_formula, ds, coef_name, args.bootstrap, rng)
        trend_rows.append(
            {
                "mutant": mutant,
                "n_observations": len(ds),
                "n_ocr": ds.index.size // 3,
                "delta_rt_beta_change_per_stage": float(fit.params[coef_name]),
                "cluster_ci_low": float(ci.iloc[0]),
                "cluster_ci_high": float(ci.iloc[1]),
                "cluster_p": float(fit.pvalues[coef_name]),
                "block_bootstrap_beta": point,
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
            }
        )
    trend = pd.DataFrame(trend_rows)
    trend.to_csv(outdir / "rt_atac_es_to_ls_trend_tests_004.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
    for j, mutant in enumerate(MUTANTS):
        a = assoc[assoc["mutant"] == mutant].set_index("stage").loc[list(STAGES)].reset_index()
        ax = axes[0, j]
        x = np.arange(len(a))
        ax.errorbar(
            x, a["standardized_beta_delta_rt"],
            yerr=[a["standardized_beta_delta_rt"] - a["bootstrap_ci_low"],
                  a["bootstrap_ci_high"] - a["standardized_beta_delta_rt"]],
            fmt="o-", color=COLORS[mutant], capsize=3,
        )
        ax.axhline(0, color="black", lw=0.7, ls="--")
        ax.set_xticks(x, STAGES)
        ax.set_ylabel("Standardized β for ΔRT")
        ax.set_title(f"{chr(97+j)}  {mutant.upper()}: adjusted ΔRT–ΔATAC", loc="left", fontweight="bold")

        ax2 = axes[1, j]
        bd = binned_df[binned_df["mutant"] == mutant]
        for stage, color in zip(STAGES, ["#4477AA", "#66AA55", "#CC6677"]):
            q = bd[bd["stage"] == stage].sort_values("mean_delta_rt")
            ax2.plot(q["mean_delta_rt"], q["mean_delta_atac"], marker="o", ms=3, label=stage, color=color)
            ax2.fill_between(q["mean_delta_rt"], q["ci_low"], q["ci_high"], color=color, alpha=0.13)
        ax2.axhline(0, color="black", lw=0.6, ls="--")
        ax2.axvline(0, color="black", lw=0.6, ls="--")
        ax2.set_xlabel("Mean ΔRT within decile")
        ax2.set_ylabel("Mean ATAC log2 fold change")
        ax2.set_title(f"{chr(99+j)}  Block-robust decile summaries", loc="left", fontweight="bold")
        ax2.legend(frameon=False, ncol=3)
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_atac_stage_association_004.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_atac_stage_association_004.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Mapped OCRs:", len(metrics))
    print(assoc.to_string(index=False))
    print("\nES-to-LS trend tests")
    print(trend.to_string(index=False))


if __name__ == "__main__":
    main()
