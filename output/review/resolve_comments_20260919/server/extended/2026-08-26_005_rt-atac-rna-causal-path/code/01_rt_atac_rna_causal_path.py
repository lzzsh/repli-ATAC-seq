#!/usr/bin/env python3
"""Estimate genetic-perturbation-anchored RT -> promoter ATAC -> RNA paths.

This is a locus-level path analysis, not proof of individual-level mediation:
Repli-seq, phased ATAC-seq, and RNA-seq were measured in different sample sets.
Inference uses 1-Mb genomic-block clustered standard errors and block bootstrap.
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
from statsmodels.stats.multitest import multipletests


SEED = 20260826
STAGES = ("ES", "MS", "LS")
MUTANTS = ("cpp8", "cpp11")
COLORS = {"cpp8": "#C25759", "cpp11": "#397CA3"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ocr-metrics", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--gene-rt", required=True)
    p.add_argument("--edger-cpp8", required=True)
    p.add_argument("--edger-cpp11", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--windows", default="1000,3000")
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def parse_gene_annotation(path: str) -> pd.DataFrame:
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
            for oi in original[left:right][keep]:
                records.append((row.gene_id, int(oi)))
    return pd.DataFrame(records, columns=["gene_id", "ocr_index"]).drop_duplicates()


def aggregate_promoter_atac(
    mapping: pd.DataFrame,
    ocr: pd.DataFrame,
    genes: pd.DataFrame,
    window: int,
) -> pd.DataFrame:
    x = mapping.merge(ocr.reset_index(names="ocr_index"), on="ocr_index", how="inner")
    count_cols = [f"{genotype}_{stage.lower()}_atac" for stage in STAGES for genotype in ("wt", *MUTANTS)]
    agg = x.groupby("gene_id", sort=False).agg(
        n_promoter_ocr=("ocr_index", "nunique"),
        **{col: (col, "sum") for col in count_cols},
    ).reset_index()
    agg["promoter_window_bp"] = window
    pseudo = 0.5
    for stage in STAGES:
        sl = stage.lower()
        agg[f"wt_{sl}_log_atac"] = np.log2(agg[f"wt_{sl}_atac"] + pseudo)
        for mutant in MUTANTS:
            agg[f"delta_atac_{mutant}_{sl}"] = np.log2(agg[f"{mutant}_{sl}_atac"] + pseudo) - np.log2(
                agg[f"wt_{sl}_atac"] + pseudo
            )
    return genes.merge(agg, on="gene_id", how="inner")


def model_sufficient(formula: str, data: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    aligned_blocks = data.loc[x.index, "block"].astype(str).to_numpy()
    block_to_index = {b: i for i, b in enumerate(blocks)}
    p = x.shape[1]
    xtx = np.zeros((len(blocks), p, p), dtype=float)
    xty = np.zeros((len(blocks), p), dtype=float)
    xa = x.to_numpy(float)
    ya = y.iloc[:, 0].to_numpy(float)
    for block in np.unique(aligned_blocks):
        take = aligned_blocks == block
        bi = block_to_index[block]
        xg = xa[take]
        yg = ya[take]
        xtx[bi] = xg.T @ xg
        xty[bi] = xg.T @ yg
    return {"xtx": xtx, "xty": xty, "columns": list(x.columns)}


def solve_sample(suff: dict[str, object], sampled: np.ndarray) -> np.ndarray:
    xtx = suff["xtx"]
    xty = suff["xty"]
    return np.linalg.pinv(xtx[sampled].sum(axis=0)) @ xty[sampled].sum(axis=0)


def bootstrap_p(values: np.ndarray) -> float:
    n = np.isfinite(values).sum()
    if n == 0:
        return np.nan
    x = values[np.isfinite(values)]
    lower = (np.sum(x <= 0) + 1) / (len(x) + 1)
    upper = (np.sum(x >= 0) + 1) / (len(x) + 1)
    return float(min(1.0, 2 * min(lower, upper)))


def fit_path(d: pd.DataFrame, nboot: int, rng: np.random.Generator) -> dict[str, float]:
    a_formula = (
        "z_delta_atac ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + "
        "z_log_n_ocr + C(chr)"
    )
    b_formula = (
        "z_rna_logfc ~ z_delta_atac + z_delta_rt + z_wt_wrt + "
        "z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)"
    )
    c_formula = (
        "z_rna_logfc ~ z_delta_rt + z_wt_wrt + z_wt_log_atac + "
        "z_baseline_rna + z_log_n_ocr + C(chr)"
    )
    fit_a = smf.ols(a_formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
    fit_b = smf.ols(b_formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
    fit_c = smf.ols(c_formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})

    a = float(fit_a.params["z_delta_rt"])
    b = float(fit_b.params["z_delta_atac"])
    direct = float(fit_b.params["z_delta_rt"])
    total = float(fit_c.params["z_delta_rt"])
    blocks = np.sort(d["block"].astype(str).unique())
    sa = model_sufficient(a_formula, d, blocks)
    sb = model_sufficient(b_formula, d, blocks)
    sc = model_sufficient(c_formula, d, blocks)
    ia = sa["columns"].index("z_delta_rt")
    ib = sb["columns"].index("z_delta_atac")
    idirect = sb["columns"].index("z_delta_rt")
    itotal = sc["columns"].index("z_delta_rt")
    boot_a = np.empty(nboot)
    boot_b = np.empty(nboot)
    boot_direct = np.empty(nboot)
    boot_total = np.empty(nboot)
    ng = len(blocks)
    for i in range(nboot):
        sampled = rng.integers(0, ng, ng)
        ba = solve_sample(sa, sampled)
        bb = solve_sample(sb, sampled)
        bc = solve_sample(sc, sampled)
        boot_a[i] = ba[ia]
        boot_b[i] = bb[ib]
        boot_direct[i] = bb[idirect]
        boot_total[i] = bc[itotal]
    boot_indirect = boot_a * boot_b

    def ci(v: np.ndarray) -> tuple[float, float]:
        lo, hi = np.nanpercentile(v, [2.5, 97.5])
        return float(lo), float(hi)

    a_lo, a_hi = ci(boot_a)
    b_lo, b_hi = ci(boot_b)
    ind_lo, ind_hi = ci(boot_indirect)
    direct_lo, direct_hi = ci(boot_direct)
    total_lo, total_hi = ci(boot_total)
    return {
        "n_genes": len(d),
        "n_blocks": d["block"].nunique(),
        "path_a_delta_rt_to_delta_atac": a,
        "path_a_cluster_p": float(fit_a.pvalues["z_delta_rt"]),
        "path_a_boot_ci_low": a_lo,
        "path_a_boot_ci_high": a_hi,
        "path_b_delta_atac_to_rna": b,
        "path_b_cluster_p": float(fit_b.pvalues["z_delta_atac"]),
        "path_b_boot_ci_low": b_lo,
        "path_b_boot_ci_high": b_hi,
        "indirect_effect_a_times_b": a * b,
        "indirect_boot_ci_low": ind_lo,
        "indirect_boot_ci_high": ind_hi,
        "indirect_boot_p": bootstrap_p(boot_indirect),
        "direct_effect_c_prime": direct,
        "direct_cluster_p": float(fit_b.pvalues["z_delta_rt"]),
        "direct_boot_ci_low": direct_lo,
        "direct_boot_ci_high": direct_hi,
        "total_effect_c": total,
        "total_cluster_p": float(fit_c.pvalues["z_delta_rt"]),
        "total_boot_ci_low": total_lo,
        "total_boot_ci_high": total_hi,
        "outcome_model_adjusted_r2": float(fit_b.rsquared_adj),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    windows = [int(x) for x in args.windows.split(",")]

    genes = parse_gene_annotation(args.gff)
    genes.to_csv(outdir / "gene_tss_from_all_diy_gff3_005.tsv", sep="\t", index=False)
    ocr = pd.read_csv(args.ocr_metrics, sep="\t")
    ocr = ocr.reset_index(drop=True)

    promoter_frames = []
    mapping_rows = []
    for window in windows:
        mapping = promoter_ocr_map(genes, ocr, window)
        mapping_rows.append({"promoter_window_bp": window, "mapped_genes": mapping["gene_id"].nunique(), "gene_ocr_pairs": len(mapping)})
        promoter_frames.append(aggregate_promoter_atac(mapping, ocr, genes, window))
    promoter = pd.concat(promoter_frames, ignore_index=True)
    pd.DataFrame(mapping_rows).to_csv(outdir / "promoter_ocr_mapping_summary_005.tsv", sep="\t", index=False)

    rt = pd.read_csv(args.gene_rt, sep="\t")
    edger = {
        "cpp8": pd.read_csv(args.edger_cpp8, index_col=0).reset_index(names="gene_id"),
        "cpp11": pd.read_csv(args.edger_cpp11, index_col=0).reset_index(names="gene_id"),
    }
    all_metrics = []
    path_rows = []
    for window in windows:
        base = promoter[promoter["promoter_window_bp"] == window].merge(rt, on="gene_id", how="inner")
        for mutant in MUTANTS:
            d0 = base.merge(edger[mutant][["gene_id", "logFC", "logCPM", "PValue", "FDR"]], on="gene_id", how="inner")
            d0["mutant"] = mutant
            d0["rna_logfc"] = pd.to_numeric(d0["logFC"], errors="coerce")
            d0["baseline_rna"] = np.log2(pd.to_numeric(d0["WT_mean"], errors="coerce") + 1.0)
            d0["delta_rt"] = pd.to_numeric(d0[f"deltaRT_{mutant}"], errors="coerce")
            d0["wt_wrt"] = pd.to_numeric(d0["WT_WRT"], errors="coerce")
            d0["block"] = d0["chr"].astype(str) + ":" + (d0["tss"] // 1_000_000).astype(str)
            d0["log_n_ocr"] = np.log1p(d0["n_promoter_ocr"])
            for stage in STAGES:
                sl = stage.lower()
                d = d0.copy()
                d["stage"] = stage
                d["delta_atac"] = pd.to_numeric(d[f"delta_atac_{mutant}_{sl}"], errors="coerce")
                d["wt_log_atac"] = pd.to_numeric(d[f"wt_{sl}_log_atac"], errors="coerce")
                needed = ["rna_logfc", "baseline_rna", "delta_rt", "wt_wrt", "delta_atac", "wt_log_atac", "log_n_ocr"]
                d = d.dropna(subset=needed).copy()
                for col in needed:
                    d[f"z_{col}"] = zscore(d[col])
                result = fit_path(d, args.bootstrap, rng)
                result.update({"mutant": mutant, "stage": stage, "promoter_window_bp": window})
                path_rows.append(result)
                keep_cols = [
                    "gene_id", "chr", "tss", "strand", "block", "mutant", "stage", "promoter_window_bp",
                    "n_promoter_ocr", "delta_rt", "wt_wrt", "delta_atac", "wt_log_atac",
                    "rna_logfc", "baseline_rna", "logCPM", "PValue", "FDR",
                ]
                all_metrics.append(d[keep_cols])

    paths = pd.DataFrame(path_rows)
    paths["indirect_fdr_across_12_tests"] = multipletests(paths["indirect_boot_p"], method="fdr_bh")[1]
    paths["indirect_and_total_same_sign"] = np.sign(paths["indirect_effect_a_times_b"]) == np.sign(paths["total_effect_c"])
    paths.to_csv(outdir / "rt_atac_rna_path_effects_005.tsv", sep="\t", index=False)
    pd.concat(all_metrics, ignore_index=True).to_csv(outdir / "gene_level_path_metrics_005.tsv.gz", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharex=True)
    for row_i, mutant in enumerate(MUTANTS):
        for col_i, window in enumerate(windows):
            ax = axes[row_i, col_i]
            q = paths[(paths["mutant"] == mutant) & (paths["promoter_window_bp"] == window)].set_index("stage").loc[list(STAGES)].reset_index()
            y = np.arange(len(q))
            ax.errorbar(
                q["indirect_effect_a_times_b"], y,
                xerr=[q["indirect_effect_a_times_b"] - q["indirect_boot_ci_low"],
                      q["indirect_boot_ci_high"] - q["indirect_effect_a_times_b"]],
                fmt="o", color=COLORS[mutant], capsize=3, label="Indirect a×b",
            )
            ax.scatter(q["direct_effect_c_prime"], y + 0.12, marker="x", color="#555555", label="Direct c′")
            ax.axvline(0, color="black", lw=0.7, ls="--")
            ax.set_yticks(y, STAGES)
            ax.set_title(
                f"{chr(97 + row_i * 2 + col_i)}  {mutant.upper()}, TSS ±{window // 1000} kb",
                loc="left", fontweight="bold",
            )
            ax.set_xlabel("Standardized path effect")
            if row_i == 0 and col_i == 0:
                ax.legend(frameon=False, loc="best")
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "rt_atac_rna_causal_paths_005.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "rt_atac_rna_causal_paths_005.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Annotation genes:", len(genes))
    print(pd.DataFrame(mapping_rows).to_string(index=False))
    print("\nPath results")
    print(paths.to_string(index=False))


if __name__ == "__main__":
    main()
