#!/usr/bin/env python3
"""Replicate the CPP8 RT-ATAC-RNA statistical path in two independent alleles.

All locus-level inference is spatially resampled by 1-Mb genomic blocks.  The
analysis tests a perturbation-anchored path and does not identify individual-
level mediation because the three assays were measured in separate samples.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrices
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests


SEED = 20260826
STAGES = ("ES", "MS", "LS")
ALLELES = ("cpp8_1", "cpp8_3")
COLORS = {"cpp8_1": "#8F3B46", "cpp8_3": "#D67B5B"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--atac-counts", required=True)
    p.add_argument("--repli", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--rna", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap-ocr", type=int, default=1000)
    p.add_argument("--bootstrap-path", type=int, default=2000)
    p.add_argument("--windows", default="1000,3000")
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def wrt(es: pd.Series, ms: pd.Series, ls: pd.Series) -> np.ndarray:
    esv = pd.to_numeric(es, errors="coerce").to_numpy(float)
    msv = pd.to_numeric(ms, errors="coerce").to_numpy(float)
    lsv = pd.to_numeric(ls, errors="coerce").to_numpy(float)
    denom = esv + msv + lsv
    return np.divide(0.5 * msv + lsv, denom, out=np.full_like(denom, np.nan), where=denom > 0)


def positive_count_size_factors(counts: np.ndarray) -> np.ndarray:
    positive = counts > 0
    npos = positive.sum(axis=1)
    logs = np.zeros_like(counts, dtype=float)
    np.log(counts, out=logs, where=positive)
    geo = np.exp(logs.sum(axis=1) / np.maximum(npos, 1))
    usable = (npos >= 2) & np.isfinite(geo) & (geo > 0)
    sf = np.empty(counts.shape[1], dtype=float)
    for j in range(counts.shape[1]):
        valid = usable & positive[:, j]
        sf[j] = np.median(counts[valid, j] / geo[valid])
    if not np.all(np.isfinite(sf) & (sf > 0)):
        raise ValueError("Non-positive ATAC size factor")
    return sf / np.exp(np.mean(np.log(sf)))


def read_rt(path: str) -> pd.DataFrame:
    names = [
        "chr", "start", "end",
        "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
        "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
        "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
        "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
        "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS",
    ]
    raw = pd.read_csv(path, sep="\t", header=None, names=names)
    raw = raw[raw["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy()
    raw[["start", "end"]] = raw[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    raw = raw.dropna(subset=["start", "end"]).reset_index(drop=True)
    raw[["start", "end"]] = raw[["start", "end"]].astype(int)
    library_cols = names[3:]
    raw[library_cols] = raw[library_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    length_kb = (raw["end"] - raw["start"] + 1).to_numpy(float) / 1000.0
    rpk = raw[library_cols].to_numpy(float) / length_kb[:, None]
    denom = rpk.sum(axis=0, keepdims=True) / 1_000_000.0
    raw[library_cols] = np.divide(rpk, denom, out=np.zeros_like(rpk), where=denom > 0)

    replicate_map = {
        "WT": {"G1": ["WT-1-G1", "WT-2-G1"], "ES": ["WT-1-ES", "WT-2-ES"], "MS": ["WT-1-MS"], "LS": ["WT-1-LS"]},
        "cpp8_1": {"G1": ["tcx2_1-1-G1", "tcx2_1-2-G1"], "ES": ["tcx2_1-1-ES", "tcx2_1-2-ES"], "MS": ["tcx2_1-1-MS"], "LS": ["tcx2_1-1-LS"]},
        "cpp8_3": {"G1": ["tcx2_3-1-G1"], "ES": ["tcx2_3-1-ES"], "MS": ["tcx2_3-1-MS"], "LS": ["tcx2_3-1-LS"]},
    }
    rt = raw[["chr", "start", "end"]].copy()
    for genotype, stage_map in replicate_map.items():
        g1 = raw[stage_map["G1"]].mean(axis=1)
        for stage in STAGES:
            signal = raw[stage_map[stage]].mean(axis=1)
            rt[f"{genotype}_{stage}_norm"] = signal / (g1 + 1e-6)
        rt[f"{genotype}_wrt"] = wrt(
            rt[f"{genotype}_ES_norm"], rt[f"{genotype}_MS_norm"], rt[f"{genotype}_LS_norm"]
        )
    for allele in ALLELES:
        rt[f"delta_rt_{allele}"] = rt[f"{allele}_wrt"] - rt["WT_wrt"]
    rt["block"] = rt["chr"].astype(str) + ":" + (rt["start"] // 1_000_000).astype(str)
    return rt


def map_intervals_to_rt(queries: pd.DataFrame, rt: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    result = np.full((len(queries), len(value_cols)), np.nan, dtype=float)
    overlap_bp = np.zeros(len(queries), dtype=int)
    n_bins = np.zeros(len(queries), dtype=int)
    for chrom, qi in queries.groupby("chr", sort=False).groups.items():
        qd = queries.loc[qi].sort_values("start")
        rd = rt[rt["chr"] == chrom].sort_values("start")
        if rd.empty:
            continue
        rs = rd["start"].to_numpy(int)
        re_ = rd["end"].to_numpy(int)
        rv = rd[value_cols].to_numpy(float)
        for idx, start, end in qd[["start", "end"]].itertuples():
            left = int(np.searchsorted(re_, start, side="right"))
            right = int(np.searchsorted(rs, end, side="left"))
            if right <= left:
                continue
            ov = np.minimum(re_[left:right], end) - np.maximum(rs[left:right], start)
            keep = ov > 0
            if not np.any(keep):
                continue
            weights = ov[keep].astype(float)
            vals = rv[left:right][keep]
            valid = np.isfinite(vals)
            den = (valid * weights[:, None]).sum(axis=0)
            num = (np.where(valid, vals, 0.0) * weights[:, None]).sum(axis=0)
            result[idx] = np.divide(num, den, out=np.full(len(value_cols), np.nan), where=den > 0)
            overlap_bp[idx] = int(weights.sum())
            n_bins[idx] = int(keep.sum())
    out = pd.DataFrame(result, columns=value_cols, index=queries.index)
    out["rt_overlap_bp"] = overlap_bp
    out["n_rt_bins"] = n_bins
    return out


def parse_genes(path: str) -> pd.DataFrame:
    names = ["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    gff = pd.read_csv(path, sep="\t", comment="#", header=None, names=names, low_memory=False)
    genes = gff[
        (gff["feature"] == "gene")
        & gff["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")
        & gff["attributes"].astype(str).str.contains(r"(?:^|;)ID=LOC_Os", regex=True)
    ].copy()
    genes["gene_id"] = genes["attributes"].astype(str).str.extract(r"(?:^|;)ID=([^;]+)", expand=False)
    genes[["start", "end"]] = genes[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    genes = genes.dropna(subset=["gene_id", "start", "end", "strand"])
    genes[["start", "end"]] = genes[["start", "end"]].astype(int)
    genes["tss"] = np.where(genes["strand"] == "+", genes["start"] - 1, genes["end"] - 1)
    genes = genes.sort_values(["gene_id", "source"]).drop_duplicates("gene_id", keep="first")
    return genes[["gene_id", "chr", "start", "end", "strand", "tss", "source"]].reset_index(drop=True)


def map_tss_to_rt(genes: pd.DataFrame, rt: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    out = genes.copy()
    for col in value_cols:
        out[col] = np.nan
    for chrom, gi in out.groupby("chr", sort=False).groups.items():
        rd = rt[rt["chr"] == chrom].sort_values("start")
        if rd.empty:
            continue
        starts = rd["start"].to_numpy(int)
        ends = rd["end"].to_numpy(int)
        vals = rd[value_cols].to_numpy(float)
        tss = out.loc[gi, "tss"].to_numpy(int)
        pos = np.searchsorted(starts, tss, side="right") - 1
        valid = (pos >= 0) & (tss < ends[np.maximum(pos, 0)])
        mapped = np.full((len(tss), len(value_cols)), np.nan)
        mapped[valid] = vals[pos[valid]]
        out.loc[gi, value_cols] = mapped
    out["block"] = out["chr"].astype(str) + ":" + (out["tss"] // 1_000_000).astype(str)
    return out


def promoter_ocr_map(genes: pd.DataFrame, ocr: pd.DataFrame, window: int) -> pd.DataFrame:
    records: list[tuple[str, int]] = []
    for chrom, gd in genes.groupby("chr", sort=False):
        od = ocr[ocr["chr"] == chrom].sort_values("start")
        if od.empty:
            continue
        starts = od["start"].to_numpy(int)
        ends = od["end"].to_numpy(int)
        original = od.index.to_numpy(int)
        max_len = int(np.max(ends - starts))
        for row in gd[["gene_id", "tss"]].itertuples(index=False):
            pstart = max(0, int(row.tss) - window)
            pend = int(row.tss) + window + 1
            left = int(np.searchsorted(starts, pstart - max_len, side="left"))
            right = int(np.searchsorted(starts, pend, side="left"))
            if right <= left:
                continue
            keep = ends[left:right] > pstart
            records.extend((row.gene_id, int(oi)) for oi in original[left:right][keep])
    return pd.DataFrame(records, columns=["gene_id", "ocr_index"]).drop_duplicates()


def model_sufficient(formula: str, data: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    aligned = data.loc[x.index, "block"].astype(str).to_numpy()
    lookup = {b: i for i, b in enumerate(blocks)}
    p = x.shape[1]
    xtx = np.zeros((len(blocks), p, p), dtype=float)
    xty = np.zeros((len(blocks), p), dtype=float)
    xa = x.to_numpy(float)
    ya = y.iloc[:, 0].to_numpy(float)
    for block in np.unique(aligned):
        take = aligned == block
        bi = lookup[block]
        xtx[bi] = xa[take].T @ xa[take]
        xty[bi] = xa[take].T @ ya[take]
    return {"xtx": xtx, "xty": xty, "columns": list(x.columns)}


def solve_sufficient(suff: dict[str, object], sampled: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(suff["xtx"][sampled].sum(axis=0)) @ suff["xty"][sampled].sum(axis=0)


def bootstrap_p(values: np.ndarray) -> float:
    x = values[np.isfinite(values)]
    if len(x) == 0:
        return np.nan
    lower = (np.sum(x <= 0) + 1) / (len(x) + 1)
    upper = (np.sum(x >= 0) + 1) / (len(x) + 1)
    return float(min(1.0, 2 * min(lower, upper)))


def block_bootstrap_coefficient(formula: str, data: pd.DataFrame, coefficient: str, nboot: int, rng: np.random.Generator) -> tuple[float, float, float]:
    blocks = np.sort(data["block"].astype(str).unique())
    suff = model_sufficient(formula, data, blocks)
    idx = suff["columns"].index(coefficient)
    sampled_all = np.arange(len(blocks))
    point = float(solve_sufficient(suff, sampled_all)[idx])
    boots = np.empty(nboot)
    for i in range(nboot):
        boots[i] = solve_sufficient(suff, rng.integers(0, len(blocks), len(blocks)))[idx]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def prepare_path(base: pd.DataFrame, allele: str, stage: str) -> pd.DataFrame:
    sl = stage.lower()
    d = base.copy()
    d["allele"] = allele
    d["stage"] = stage
    d["delta_rt"] = pd.to_numeric(d[f"delta_rt_{allele}"], errors="coerce")
    d["delta_atac"] = pd.to_numeric(d[f"delta_atac_{allele}_{sl}"], errors="coerce")
    d["rna_logfc"] = pd.to_numeric(d[f"{allele}_logFC"], errors="coerce")
    d["wt_wrt"] = pd.to_numeric(d["WT_wrt"], errors="coerce")
    d["wt_log_atac"] = pd.to_numeric(d[f"wt_{sl}_log_atac"], errors="coerce")
    d["baseline_rna"] = np.log2(pd.to_numeric(d["wt_mean_cpm"], errors="coerce") + 1.0)
    d["log_n_ocr"] = np.log1p(pd.to_numeric(d["n_promoter_ocr"], errors="coerce"))
    needed = ["delta_rt", "delta_atac", "rna_logfc", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]
    d = d.dropna(subset=needed).copy()
    for col in needed:
        d[f"z_{col}"] = zscore(d[col])
    return d


PATH_FORMULAS = {
    "a": "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + z_log_n_ocr + C(chr)",
    "b": "z_rna_logfc ~ z_delta_atac + z_delta_rt + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
    "c": "z_rna_logfc ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)",
}


def path_point_and_sufficient(d: pd.DataFrame, blocks: np.ndarray) -> tuple[dict[str, float], dict[str, dict[str, object]]]:
    fits = {
        key: smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
        for key, formula in PATH_FORMULAS.items()
    }
    suff = {key: model_sufficient(formula, d, blocks) for key, formula in PATH_FORMULAS.items()}
    a = float(fits["a"].params["z_delta_rt"])
    b = float(fits["b"].params["z_delta_atac"])
    direct = float(fits["b"].params["z_delta_rt"])
    total = float(fits["c"].params["z_delta_rt"])
    point = {
        "a": a, "b": b, "indirect": a * b, "direct": direct, "total": total,
        "a_cluster_p": float(fits["a"].pvalues["z_delta_rt"]),
        "b_cluster_p": float(fits["b"].pvalues["z_delta_atac"]),
        "direct_cluster_p": float(fits["b"].pvalues["z_delta_rt"]),
        "total_cluster_p": float(fits["c"].pvalues["z_delta_rt"]),
        "outcome_adjusted_r2": float(fits["b"].rsquared_adj),
    }
    return point, suff


def extract_path_boot(suff: dict[str, dict[str, object]], sampled: np.ndarray) -> dict[str, float]:
    ba = solve_sufficient(suff["a"], sampled)
    bb = solve_sufficient(suff["b"], sampled)
    bc = solve_sufficient(suff["c"], sampled)
    a = float(ba[suff["a"]["columns"].index("z_delta_rt")])
    b = float(bb[suff["b"]["columns"].index("z_delta_atac")])
    return {
        "a": a,
        "b": b,
        "indirect": a * b,
        "direct": float(bb[suff["b"]["columns"].index("z_delta_rt")]),
        "total": float(bc[suff["c"]["columns"].index("z_delta_rt")]),
    }


def ci(values: np.ndarray) -> tuple[float, float]:
    lo, hi = np.nanpercentile(values, [2.5, 97.5])
    return float(lo), float(hi)


def block_bootstrap_correlation(d: pd.DataFrame, xcol: str, ycol: str, nboot: int, rng: np.random.Generator) -> dict[str, float]:
    q = d[[xcol, ycol, "block"]].dropna().reset_index(drop=True)
    blocks = np.array(sorted(q["block"].astype(str).unique()))
    indices = {b: np.flatnonzero(q["block"].astype(str).to_numpy() == b) for b in blocks}
    x = q[xcol].to_numpy(float)
    y = q[ycol].to_numpy(float)
    pear = float(pearsonr(x, y).statistic)
    spear = float(spearmanr(x, y).statistic)
    bp = np.empty(nboot)
    bs = np.empty(nboot)
    for i in range(nboot):
        sampled = rng.integers(0, len(blocks), len(blocks))
        take = np.concatenate([indices[blocks[j]] for j in sampled])
        bp[i] = pearsonr(x[take], y[take]).statistic
        bs[i] = spearmanr(x[take], y[take]).statistic
    pl, ph = ci(bp)
    sl, sh = ci(bs)
    return {"n": len(q), "n_blocks": len(blocks), "pearson_r": pear, "pearson_ci_low": pl, "pearson_ci_high": ph, "spearman_rho": spear, "spearman_ci_low": sl, "spearman_ci_high": sh}


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    windows = [int(x) for x in args.windows.split(",")]

    atac = pd.read_csv(args.atac_counts, sep="\t")
    atac = atac[atac["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy().reset_index(drop=True)
    atac[["start", "end"]] = atac[["start", "end"]].astype(int)
    sample_cols = [c for c in atac.columns if c not in {"chr", "start", "end"}]
    raw_counts = atac[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float)
    size_factors = positive_count_size_factors(raw_counts)
    atac[sample_cols] = raw_counts / size_factors[None, :]
    pd.DataFrame({"sample": sample_cols, "size_factor": size_factors}).to_csv(outdir / "atac_positive_count_size_factors_009.tsv", sep="\t", index=False)

    group_cols = {
        "ES": {"wt": ["WT_ES_1", "WT_ES_2"], "cpp8_1": ["CPP8_1_ES_1", "CPP8_1_ES_2"], "cpp8_3": ["CPP8_3_ES_1", "CPP8_3_ES_2"]},
        "MS": {"wt": ["WT_MS"], "cpp8_1": ["CPP8_1_MS"], "cpp8_3": ["CPP8_3_MS"]},
        "LS": {"wt": ["WT_LS"], "cpp8_1": ["CPP8_1_LS"], "cpp8_3": ["CPP8_3_LS"]},
    }
    for stage, groups in group_cols.items():
        sl = stage.lower()
        for genotype, cols in groups.items():
            atac[f"{genotype}_{sl}_atac"] = atac[cols].mean(axis=1)
        atac[f"wt_{sl}_log_atac"] = np.log2(atac[f"wt_{sl}_atac"] + 0.5)
        for allele in ALLELES:
            atac[f"delta_atac_{allele}_{sl}"] = np.log2(atac[f"{allele}_{sl}_atac"] + 0.5) - np.log2(atac[f"wt_{sl}_atac"] + 0.5)

    rt = read_rt(args.repli)
    rt.to_csv(outdir / "repli_rt_bin_metrics_009.tsv.gz", sep="\t", index=False)
    value_cols = [
        "WT_wrt", "WT_ES_norm", "WT_MS_norm", "WT_LS_norm",
        "delta_rt_cpp8_1", "delta_rt_cpp8_3",
    ]
    mapped = map_intervals_to_rt(atac[["chr", "start", "end"]], rt, value_cols)
    ocr = pd.concat([atac, mapped], axis=1)
    ocr["length"] = ocr["end"] - ocr["start"]
    ocr["block"] = ocr["chr"].astype(str) + ":" + (ocr["start"] // 1_000_000).astype(str)
    ocr = ocr[(ocr["rt_overlap_bp"] > 0) & ocr["WT_wrt"].notna()].reset_index(drop=True)

    ocr_rows = []
    ocr_formula = "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + z_log_length + C(chr)"
    for allele in ALLELES:
        for stage in STAGES:
            sl = stage.lower()
            d = ocr[["chr", "block", "length", "WT_wrt", f"delta_rt_{allele}", f"wt_{sl}_log_atac", f"delta_atac_{allele}_{sl}"]].copy()
            d.columns = ["chr", "block", "length", "wt_wrt", "delta_rt", "wt_log_atac", "delta_atac"]
            d = d.dropna().copy()
            d["log_length"] = np.log1p(d["length"])
            for col in ["delta_atac", "delta_rt", "wt_wrt", "wt_log_atac", "log_length"]:
                d[f"z_{col}"] = zscore(d[col])
            fit = smf.ols(ocr_formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
            point, lo, hi = block_bootstrap_coefficient(ocr_formula, d, "z_delta_rt", args.bootstrap_ocr, rng)
            ocr_rows.append({
                "allele": allele, "stage": stage, "n_ocr": len(d), "n_blocks": d["block"].nunique(),
                "standardized_beta_delta_rt": float(fit.params["z_delta_rt"]),
                "cluster_se": float(fit.bse["z_delta_rt"]), "cluster_p": float(fit.pvalues["z_delta_rt"]),
                "block_bootstrap_beta": point, "bootstrap_ci_low": lo, "bootstrap_ci_high": hi,
                "adjusted_r2": float(fit.rsquared_adj),
            })
    ocr_assoc = pd.DataFrame(ocr_rows)
    ocr_assoc["fdr_across_six_tests"] = multipletests(ocr_assoc["cluster_p"], method="fdr_bh")[1]
    ocr_assoc.to_csv(outdir / "cpp8_allele_ocr_associations_009.tsv", sep="\t", index=False)
    ocr.to_csv(outdir / "ocr_rt_atac_metrics_009.tsv.gz", sep="\t", index=False)

    genes = parse_genes(args.gff)
    gene_rt = map_tss_to_rt(genes, rt, value_cols)
    gene_rt.to_csv(outdir / "gene_tss_rt_metrics_009.tsv.gz", sep="\t", index=False)
    rna = pd.read_csv(args.rna, sep="\t")

    promoter_frames = []
    map_rows = []
    count_cols = [f"{g}_{s.lower()}_atac" for s in STAGES for g in ("wt", *ALLELES)]
    for window in windows:
        mapping = promoter_ocr_map(genes, ocr, window)
        x = mapping.merge(ocr.reset_index(names="ocr_index"), on="ocr_index", how="inner")
        agg = x.groupby("gene_id", sort=False).agg(
            n_promoter_ocr=("ocr_index", "nunique"),
            **{col: (col, "sum") for col in count_cols},
        ).reset_index()
        agg["promoter_window_bp"] = window
        for stage in STAGES:
            sl = stage.lower()
            agg[f"wt_{sl}_log_atac"] = np.log2(agg[f"wt_{sl}_atac"] + 0.5)
            for allele in ALLELES:
                agg[f"delta_atac_{allele}_{sl}"] = np.log2(agg[f"{allele}_{sl}_atac"] + 0.5) - np.log2(agg[f"wt_{sl}_atac"] + 0.5)
        promoter_frames.append(agg)
        map_rows.append({"promoter_window_bp": window, "mapped_genes": agg["gene_id"].nunique(), "gene_ocr_pairs": len(mapping)})
    promoter = pd.concat(promoter_frames, ignore_index=True)
    pd.DataFrame(map_rows).to_csv(outdir / "promoter_ocr_mapping_summary_009.tsv", sep="\t", index=False)

    path_rows = []
    hetero_rows = []
    metric_frames = []
    concordance_rows = []
    for window in windows:
        base = promoter[promoter["promoter_window_bp"] == window].merge(gene_rt, on="gene_id", how="inner").merge(rna, on="gene_id", how="inner")
        for stage in STAGES:
            prepared = {allele: prepare_path(base, allele, stage) for allele in ALLELES}
            common_genes = set(prepared[ALLELES[0]]["gene_id"]) & set(prepared[ALLELES[1]]["gene_id"])
            prepared = {allele: prepared[allele][prepared[allele]["gene_id"].isin(common_genes)].sort_values("gene_id").reset_index(drop=True) for allele in ALLELES}
            if not np.array_equal(prepared[ALLELES[0]]["gene_id"].to_numpy(), prepared[ALLELES[1]]["gene_id"].to_numpy()):
                raise RuntimeError("Allele gene universes are not aligned")
            blocks = np.array(sorted(set(prepared[ALLELES[0]]["block"]) | set(prepared[ALLELES[1]]["block"])))
            points = {}
            suffs = {}
            for allele in ALLELES:
                points[allele], suffs[allele] = path_point_and_sufficient(prepared[allele], blocks)
            boot = {allele: {k: np.empty(args.bootstrap_path) for k in ["a", "b", "indirect", "direct", "total"]} for allele in ALLELES}
            for i in range(args.bootstrap_path):
                sampled = rng.integers(0, len(blocks), len(blocks))
                for allele in ALLELES:
                    vals = extract_path_boot(suffs[allele], sampled)
                    for key, value in vals.items():
                        boot[allele][key][i] = value
            for allele in ALLELES:
                row = {
                    "allele": allele, "stage": stage, "promoter_window_bp": window,
                    "n_genes": len(prepared[allele]), "n_blocks": len(blocks),
                    "path_a_delta_rt_to_delta_atac": points[allele]["a"],
                    "path_a_cluster_p": points[allele]["a_cluster_p"],
                    "path_b_delta_atac_to_rna": points[allele]["b"],
                    "path_b_cluster_p": points[allele]["b_cluster_p"],
                    "indirect_effect_a_times_b": points[allele]["indirect"],
                    "indirect_boot_p": bootstrap_p(boot[allele]["indirect"]),
                    "direct_effect_c_prime": points[allele]["direct"],
                    "direct_cluster_p": points[allele]["direct_cluster_p"],
                    "total_effect_c": points[allele]["total"],
                    "total_cluster_p": points[allele]["total_cluster_p"],
                    "outcome_model_adjusted_r2": points[allele]["outcome_adjusted_r2"],
                }
                for key, label in [("a", "path_a"), ("b", "path_b"), ("indirect", "indirect"), ("direct", "direct"), ("total", "total")]:
                    lo, hi = ci(boot[allele][key])
                    row[f"{label}_boot_ci_low"] = lo
                    row[f"{label}_boot_ci_high"] = hi
                path_rows.append(row)
                keep = ["gene_id", "chr", "tss", "strand", "block", "allele", "stage", "promoter_window_bp", "n_promoter_ocr", "delta_rt", "wt_wrt", "delta_atac", "wt_log_atac", "rna_logfc", "baseline_rna"]
                metric_frames.append(prepared[allele][keep])

            for key, label in [("a", "path_a"), ("b", "path_b"), ("indirect", "indirect"), ("direct", "direct"), ("total", "total")]:
                diff = boot["cpp8_1"][key] - boot["cpp8_3"][key]
                mean = (boot["cpp8_1"][key] + boot["cpp8_3"][key]) / 2
                dlo, dhi = ci(diff)
                mlo, mhi = ci(mean)
                hetero_rows.append({
                    "stage": stage, "promoter_window_bp": window, "effect": label,
                    "cpp8_1_point": points["cpp8_1"][key], "cpp8_3_point": points["cpp8_3"][key],
                    "difference_cpp8_1_minus_cpp8_3": points["cpp8_1"][key] - points["cpp8_3"][key],
                    "difference_boot_ci_low": dlo, "difference_boot_ci_high": dhi,
                    "difference_boot_p": bootstrap_p(diff),
                    "mean_allele_effect": (points["cpp8_1"][key] + points["cpp8_3"][key]) / 2,
                    "mean_boot_ci_low": mlo, "mean_boot_ci_high": mhi, "mean_boot_p": bootstrap_p(mean),
                })

            pair = prepared["cpp8_1"][["gene_id", "block", "delta_atac"]].rename(columns={"delta_atac": "cpp8_1"}).merge(
                prepared["cpp8_3"][["gene_id", "delta_atac"]].rename(columns={"delta_atac": "cpp8_3"}), on="gene_id"
            )
            corr = block_bootstrap_correlation(pair, "cpp8_1", "cpp8_3", args.bootstrap_ocr, rng)
            corr.update({"metric": "promoter_delta_atac", "stage": stage, "promoter_window_bp": window})
            concordance_rows.append(corr)

    paths = pd.DataFrame(path_rows)
    paths["indirect_fdr_across_12_tests"] = multipletests(paths["indirect_boot_p"], method="fdr_bh")[1]
    paths.to_csv(outdir / "cpp8_allele_path_effects_009.tsv", sep="\t", index=False)
    hetero = pd.DataFrame(hetero_rows)
    hetero["difference_fdr_within_effect"] = hetero.groupby("effect")["difference_boot_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    hetero.to_csv(outdir / "cpp8_allele_path_heterogeneity_009.tsv", sep="\t", index=False)
    metrics = pd.concat(metric_frames, ignore_index=True)
    metrics.to_csv(outdir / "gene_level_path_metrics_009.tsv.gz", sep="\t", index=False)

    rt_corr = block_bootstrap_correlation(rt, "delta_rt_cpp8_1", "delta_rt_cpp8_3", args.bootstrap_ocr, rng)
    rt_corr.update({"metric": "rt_bin_delta_rt", "stage": "all", "promoter_window_bp": 0})
    concordance_rows.append(rt_corr)
    rna_pair = metrics[(metrics["stage"] == "MS") & (metrics["promoter_window_bp"] == windows[0])]
    rna_pair = rna_pair.pivot_table(index=["gene_id", "chr", "tss", "block"], columns="allele", values="rna_logfc").reset_index()
    rna_corr = block_bootstrap_correlation(rna_pair, "cpp8_1", "cpp8_3", args.bootstrap_ocr, rng)
    rna_corr.update({"metric": "rna_logfc", "stage": "all", "promoter_window_bp": 0})
    concordance_rows.append(rna_corr)
    concordance = pd.DataFrame(concordance_rows)
    concordance.to_csv(outdir / "cpp8_allele_concordance_009.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    axes[0, 0].hexbin(rt["delta_rt_cpp8_1"], rt["delta_rt_cpp8_3"], gridsize=65, mincnt=1, cmap="magma", bins="log")
    lim = np.nanpercentile(np.abs(rt[["delta_rt_cpp8_1", "delta_rt_cpp8_3"]].to_numpy()), 99)
    axes[0, 0].plot([-lim, lim], [-lim, lim], ls="--", lw=0.8, color="white")
    axes[0, 0].set(xlabel="oscpp8-1 ΔRT", ylabel="oscpp8-3 ΔRT", title=f"a  Allelic RT concordance (r={rt_corr['pearson_r']:.3f})")

    for i, allele in enumerate(ALLELES):
        q = ocr_assoc[ocr_assoc["allele"] == allele].set_index("stage").loc[list(STAGES)].reset_index()
        y = np.arange(3) + (i - 0.5) * 0.16
        axes[0, 1].errorbar(q["standardized_beta_delta_rt"], y, xerr=[q["standardized_beta_delta_rt"] - q["bootstrap_ci_low"], q["bootstrap_ci_high"] - q["standardized_beta_delta_rt"]], fmt="o", capsize=3, color=COLORS[allele], label=allele.replace("cpp8_", "oscpp8-"))
    axes[0, 1].axvline(0, color="black", lw=0.7, ls="--")
    axes[0, 1].set_yticks(range(3), STAGES)
    axes[0, 1].set_xlabel("Standardized β: ΔRT → ΔATAC")
    axes[0, 1].set_title("b  OCR-level association")
    axes[0, 1].legend(frameon=False)

    q = paths.copy()
    q["label"] = q["stage"] + " ±" + (q["promoter_window_bp"] // 1000).astype(str) + " kb"
    order = [f"{s} ±{w // 1000} kb" for w in windows for s in STAGES]
    pos = {x: i for i, x in enumerate(order)}
    for i, allele in enumerate(ALLELES):
        qa = q[q["allele"] == allele]
        y = np.array([pos[x] for x in qa["label"]]) + (i - 0.5) * 0.16
        axes[1, 0].errorbar(qa["indirect_effect_a_times_b"], y, xerr=[qa["indirect_effect_a_times_b"] - qa["indirect_boot_ci_low"], qa["indirect_boot_ci_high"] - qa["indirect_effect_a_times_b"]], fmt="o", capsize=3, color=COLORS[allele], label=allele.replace("cpp8_", "oscpp8-"))
    axes[1, 0].axvline(0, color="black", lw=0.7, ls="--")
    axes[1, 0].set_yticks(range(len(order)), order)
    axes[1, 0].set_xlabel("Standardized indirect effect a×b")
    axes[1, 0].set_title("c  Allele-specific promoter paths")

    axes[1, 1].hexbin(rna_pair["cpp8_1"], rna_pair["cpp8_3"], gridsize=55, mincnt=1, cmap="viridis", bins="log")
    axes[1, 1].axhline(0, color="white", lw=0.5, ls="--")
    axes[1, 1].axvline(0, color="white", lw=0.5, ls="--")
    axes[1, 1].set(xlabel="oscpp8-1 RNA logFC", ylabel="oscpp8-3 RNA logFC", title=f"d  Transcriptional concordance (r={rna_corr['pearson_r']:.3f})")
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "cpp8_allele_path_replication_009.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "cpp8_allele_path_replication_009.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("OCRs mapped to RT:", len(ocr))
    print("\nOCR associations")
    print(ocr_assoc.to_string(index=False))
    print("\nPath effects")
    print(paths.to_string(index=False))
    print("\nAllele concordance")
    print(concordance.to_string(index=False))


if __name__ == "__main__":
    main()
