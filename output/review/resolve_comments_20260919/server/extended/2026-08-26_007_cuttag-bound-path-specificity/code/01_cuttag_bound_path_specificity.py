#!/usr/bin/env python3
"""Test whether RT-ATAC-RNA paths are stronger at CPP CUT&Tag-bound promoters."""

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
COLORS = {"cpp8": "#C25759", "cpp11": "#397CA3"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-metrics", required=True)
    p.add_argument("--cpp8-peaks", required=True)
    p.add_argument("--cpp11-peaks", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-dir", required=True)
    p.add_argument("--bootstrap", type=int, default=1000)
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def load_peaks(path: str) -> pd.DataFrame:
    p = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    p = p[p["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy()
    p[["start", "end"]] = p[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    return p.dropna().astype({"start": int, "end": int}).drop_duplicates()


def assign_bound(genes: pd.DataFrame, peaks: pd.DataFrame, window: int) -> pd.Series:
    bound = pd.Series(0, index=genes.index, dtype=int)
    for chrom, idx in genes.groupby("chr", sort=False).groups.items():
        p = peaks[peaks["chr"] == chrom].sort_values("start")
        if p.empty:
            continue
        ps = p["start"].to_numpy(int)
        pe = p["end"].to_numpy(int)
        max_len = int(np.max(pe - ps))
        for i in idx:
            tss = int(genes.at[i, "tss"])
            start = max(0, tss - window)
            end = tss + window + 1
            left = int(np.searchsorted(ps, start - max_len, side="left"))
            right = int(np.searchsorted(ps, end, side="left"))
            if right > left and np.any(pe[left:right] > start):
                bound.at[i] = 1
    return bound


def sufficient(formula: str, d: pd.DataFrame, blocks: np.ndarray) -> dict[str, object]:
    y, x = dmatrices(formula, d, return_type="dataframe", NA_action="drop")
    block_map = {b: i for i, b in enumerate(blocks)}
    p = x.shape[1]
    xtx = np.zeros((len(blocks), p, p), dtype=float)
    xty = np.zeros((len(blocks), p), dtype=float)
    xa = x.to_numpy(float)
    ya = y.iloc[:, 0].to_numpy(float)
    groups = d.loc[x.index, "block"].astype(str).to_numpy()
    for block in np.unique(groups):
        take = groups == block
        bi = block_map[block]
        xtx[bi] = xa[take].T @ xa[take]
        xty[bi] = xa[take].T @ ya[take]
    return {"xtx": xtx, "xty": xty, "columns": list(x.columns)}


def solve(s: dict[str, object], sampled: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(s["xtx"][sampled].sum(axis=0)) @ s["xty"][sampled].sum(axis=0)


def bootstrap_p(v: np.ndarray) -> float:
    lower = (np.sum(v <= 0) + 1) / (len(v) + 1)
    upper = (np.sum(v >= 0) + 1) / (len(v) + 1)
    return float(min(1.0, 2 * min(lower, upper)))


def fit_interaction_path(d: pd.DataFrame, nboot: int, rng: np.random.Generator) -> dict[str, float]:
    a_formula = (
        "z_delta_atac ~ z_delta_rt * bound + z_wt_wrt + z_wt_log_atac + "
        "z_log_n_ocr + C(chr)"
    )
    b_formula = (
        "z_rna_logfc ~ z_delta_atac * bound + z_delta_rt * bound + z_wt_wrt + "
        "z_wt_log_atac + z_baseline_rna + z_log_n_ocr + C(chr)"
    )
    fa = smf.ols(a_formula, d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
    fb = smf.ols(b_formula, d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
    a0 = float(fa.params["z_delta_rt"])
    a1 = a0 + float(fa.params["z_delta_rt:bound"])
    b0 = float(fb.params["z_delta_atac"])
    b1 = b0 + float(fb.params["z_delta_atac:bound"])
    indirect0 = a0 * b0
    indirect1 = a1 * b1
    contrast = indirect1 - indirect0

    blocks = np.sort(d["block"].astype(str).unique())
    sa = sufficient(a_formula, d, blocks)
    sb = sufficient(b_formula, d, blocks)
    ia0 = sa["columns"].index("z_delta_rt")
    iai = sa["columns"].index("z_delta_rt:bound")
    ib0 = sb["columns"].index("z_delta_atac")
    ibi = sb["columns"].index("z_delta_atac:bound")
    boot0 = np.empty(nboot)
    boot1 = np.empty(nboot)
    ng = len(blocks)
    for i in range(nboot):
        sampled = rng.integers(0, ng, ng)
        ba = solve(sa, sampled)
        bb = solve(sb, sampled)
        boot0[i] = ba[ia0] * bb[ib0]
        boot1[i] = (ba[ia0] + ba[iai]) * (bb[ib0] + bb[ibi])
    diff = boot1 - boot0
    lo0, hi0 = np.percentile(boot0, [2.5, 97.5])
    lo1, hi1 = np.percentile(boot1, [2.5, 97.5])
    lod, hid = np.percentile(diff, [2.5, 97.5])
    return {
        "n_genes": len(d),
        "n_bound_genes": int(d["bound"].sum()),
        "path_a_unbound": a0,
        "path_a_bound": a1,
        "path_a_bound_interaction_cluster_p": float(fa.pvalues["z_delta_rt:bound"]),
        "path_b_unbound": b0,
        "path_b_bound": b1,
        "path_b_bound_interaction_cluster_p": float(fb.pvalues["z_delta_atac:bound"]),
        "indirect_unbound": indirect0,
        "indirect_unbound_ci_low": float(lo0),
        "indirect_unbound_ci_high": float(hi0),
        "indirect_bound": indirect1,
        "indirect_bound_ci_low": float(lo1),
        "indirect_bound_ci_high": float(hi1),
        "bound_minus_unbound_indirect": contrast,
        "contrast_ci_low": float(lod),
        "contrast_ci_high": float(hid),
        "contrast_boot_p": bootstrap_p(diff),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    plotdir = Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    metrics = pd.read_csv(args.gene_metrics, sep="\t")
    metrics["log_n_ocr"] = np.log1p(metrics["n_promoter_ocr"])
    peaks = {"cpp8": load_peaks(args.cpp8_peaks), "cpp11": load_peaks(args.cpp11_peaks)}

    bound_maps = []
    for (mutant, window), d in metrics.groupby(["mutant", "promoter_window_bp"], sort=False):
        genes = d[["gene_id", "chr", "tss"]].drop_duplicates("gene_id").reset_index(drop=True)
        genes["mutant"] = mutant
        genes["promoter_window_bp"] = int(window)
        genes["bound"] = assign_bound(genes, peaks[mutant], int(window))
        bound_maps.append(genes)
    bound_map = pd.concat(bound_maps, ignore_index=True)
    bound_map.to_csv(outdir / "cuttag_bound_promoter_map_007.tsv.gz", sep="\t", index=False)

    rows = []
    for (mutant, stage, window), d in metrics.groupby(["mutant", "stage", "promoter_window_bp"], sort=False):
        bm = bound_map[(bound_map.mutant == mutant) & (bound_map.promoter_window_bp == window)][["gene_id", "bound"]]
        d = d.merge(bm, on="gene_id", how="inner").dropna().reset_index(drop=True)
        for col in ["delta_rt", "delta_atac", "rna_logfc", "wt_wrt", "wt_log_atac", "baseline_rna", "log_n_ocr"]:
            d[f"z_{col}"] = zscore(d[col])
        result = fit_interaction_path(d, args.bootstrap, rng)
        result.update({"mutant": mutant, "stage": stage, "promoter_window_bp": int(window)})
        rows.append(result)
    results = pd.DataFrame(rows)
    results["contrast_fdr_across_12_tests"] = multipletests(results["contrast_boot_p"], method="fdr_bh")[1]
    results.to_csv(outdir / "cuttag_bound_path_specificity_007.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharex=True)
    for ax, mutant in zip(axes, ["cpp8", "cpp11"]):
        q = results[results.mutant == mutant].copy()
        q["label"] = q["stage"] + " ±" + (q["promoter_window_bp"] // 1000).astype(str) + "kb"
        q = q.sort_values(["promoter_window_bp", "stage"])
        y = np.arange(len(q))
        ax.errorbar(
            q["bound_minus_unbound_indirect"], y,
            xerr=[q["bound_minus_unbound_indirect"] - q["contrast_ci_low"], q["contrast_ci_high"] - q["bound_minus_unbound_indirect"]],
            fmt="o", color=COLORS[mutant], capsize=3,
        )
        ax.axvline(0, color="black", lw=0.7, ls="--")
        ax.set_yticks(y, q["label"])
        ax.set_xlabel("Bound − unbound indirect effect")
        ax.set_title(f"{mutant.upper()} CUT&Tag promoter specificity", loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "cuttag_bound_path_specificity_007.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "cuttag_bound_path_specificity_007.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
