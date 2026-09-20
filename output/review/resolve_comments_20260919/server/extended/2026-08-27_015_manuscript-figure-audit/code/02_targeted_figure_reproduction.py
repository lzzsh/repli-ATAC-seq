#!/usr/bin/env python3
"""Targeted numerical checks for high-risk manuscript panels.

This script does not recreate the manuscript artwork. It independently checks
the numerical inputs and selection logic behind Figures 1, 2, 3, 7 and 8 from
read-only project sources and previously completed independent analyses.
"""

from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom


PROJECTS = Path("/storage2/liuxiaodongLab/liaozizhuo/Projects")
TASK = PROJECTS / "extended/2026-08-27_015_manuscript-figure-audit"
OUT = TASK / "output"


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False)


def figure1() -> None:
    source = (
        PROJECTS
        / "repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/replication_classification_results.csv"
    )
    data = pd.read_csv(source)
    final_candidates = [
        "Final_Classification",
        "final_classification",
        "Final_classification",
        "final_phase",
        "Final",
    ]
    final_col = next((x for x in final_candidates if x in data.columns), data.columns[-1])
    counts = data[final_col].astype(str).value_counts(dropna=False)
    result = pd.DataFrame(
        {
            "rt_class": counts.index,
            "n_source": counts.values,
            "fraction_source": counts.values / len(data),
            "manuscript_total": 62468,
            "source_total": len(data),
            "difference_source_minus_manuscript": len(data) - 62468,
            "source_path": str(source),
        }
    )
    write_tsv(result, "figure1_ocr_count_check_015.tsv")


def figure2() -> None:
    density_source = PROJECTS / "GWAS/rt_density_summary.tsv"
    maf_source = PROJECTS / "GWAS/rt_maf_spectrum.tsv"
    density = pd.read_csv(density_source, sep="\t")
    density["source_path"] = str(density_source)
    maf = pd.read_csv(maf_source, sep="\t")
    maf["source_path"] = str(maf_source)
    write_tsv(density, "figure2_variant_density_check_015.tsv")
    write_tsv(maf, "figureS2_maf_spectrum_check_015.tsv")


def figure3() -> None:
    independent = (
        PROJECTS
        / "extended/2026-08-26_003_mutant-rt-robustness/output/mutant_polarization_effects_003.tsv"
    )
    effects = pd.read_csv(independent, sep="\t")
    effects["independent_source_path"] = str(independent)
    effects["supports_positive_polarization"] = effects["slope_ci_low"] > 0
    write_tsv(effects, "figure3_polarization_check_015.tsv")

    old_dir = PROJECTS / "repli-seq-CR/segmentation/RT_quintile_boxplot_results"
    rows = []
    for path in sorted(old_dir.glob("*_stats.tsv")):
        table = pd.read_csv(path, sep="\t")
        q1 = float(table.loc[table["wt_quantile"] == "Q1", "mean_diff"].iloc[0])
        q5 = float(table.loc[table["wt_quantile"] == "Q5", "mean_diff"].iloc[0])
        rows.append(
            {
                "comparison": path.name.removesuffix("_RT_quintile_boxplot_all_bins_stats.tsv"),
                "q1_mean_delta": q1,
                "q5_mean_delta": q5,
                "q5_minus_q1": q5 - q1,
                "source_path": str(path),
            }
        )
    write_tsv(pd.DataFrame(rows), "figure3_legacy_quintile_check_015.tsv")


def figure7() -> None:
    base = PROJECTS / "ATAC-seq-CR-2/diffbind"
    paths = {
        "oscpp8-1": base / "TCX2_1_vs_WT_edgeR_sig.bed",
        "oscpp8-3": base / "TCX2_3_vs_WT_edgeR_sig.bed",
        "oscpp11-5": base / "sol1_5_vs_WT_edgeR_sig.bed",
        "oscpp11-8": base / "sol1_8_vs_WT_edgeR_sig.bed",
    }
    peak_sets: dict[str, set[tuple[str, int, int]]] = {}
    counts = []
    for allele, path in paths.items():
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "strand", "log2FC"],
        )
        peak_sets[allele] = set(zip(frame["chr"], frame["start"], frame["end"]))
        counts.append(
            {
                "allele": allele,
                "n_dar": len(frame),
                "n_dar_abs_log2fc_gt_0.3": int((frame["log2FC"].abs() > 0.3).sum()),
                "n_dar_abs_log2fc_le_0.3": int((frame["log2FC"].abs() <= 0.3).sum()),
                "n_up": int((frame["log2FC"] > 0).sum()),
                "n_down": int((frame["log2FC"] < 0).sum()),
                "source_path": str(path),
            }
        )
    write_tsv(pd.DataFrame(counts), "figure7_dar_counts_015.tsv")

    overlaps = []
    for a, b in combinations(peak_sets, 2):
        intersection = len(peak_sets[a] & peak_sets[b])
        union = len(peak_sets[a] | peak_sets[b])
        overlaps.append(
            {
                "set_a": a,
                "set_b": b,
                "n_a": len(peak_sets[a]),
                "n_b": len(peak_sets[b]),
                "intersection": intersection,
                "union": union,
                "jaccard": intersection / union,
                "overlap_coefficient": intersection / min(len(peak_sets[a]), len(peak_sets[b])),
            }
        )
    overlaps.append(
        {
            "set_a": "all_four",
            "set_b": "all_four",
            "n_a": math.nan,
            "n_b": math.nan,
            "intersection": len(set.intersection(*peak_sets.values())),
            "union": len(set.union(*peak_sets.values())),
            "jaccard": len(set.intersection(*peak_sets.values())) / len(set.union(*peak_sets.values())),
            "overlap_coefficient": math.nan,
        }
    )
    write_tsv(pd.DataFrame(overlaps), "figure7_peak_overlap_check_015.tsv")


def figure8() -> None:
    atac_path = (
        PROJECTS
        / "ATAC-seq-CR-2/diffbind/TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv"
    )
    rna_path = PROJECTS / "RNA-seq/star_all_rawdata/featureCounts_out/edgeR_TCX2-3-KO_vs_WT.csv"
    atac = pd.read_csv(atac_path)
    rna = pd.read_csv(rna_path).rename(columns={"Unnamed: 0": "geneId"})

    merged = atac.merge(rna[["geneId", "logFC"]], on="geneId", how="inner")
    row_r = merged[["log2FC", "logFC"]].corr().iloc[0, 1]
    by_gene = merged.groupby("geneId", as_index=False).agg(
        log2FC=("log2FC", "mean"), logFC=("logFC", "first")
    )
    gene_r = by_gene[["log2FC", "logFC"]].corr().iloc[0, 1]
    regression = pd.DataFrame(
        [
            {
                "analysis": "manuscript_row_level",
                "n_records": len(merged),
                "n_unique_genes": merged["geneId"].nunique(),
                "n_genes_with_multiple_records": int((merged["geneId"].value_counts() > 1).sum()),
                "pearson_r": row_r,
                "r_squared": row_r**2,
                "selection": "FDR-significant RNA genes intersected with CPP8-bound promoter-proximal significant DAR records",
            },
            {
                "analysis": "mean_atac_per_unique_gene_sensitivity",
                "n_records": len(by_gene),
                "n_unique_genes": len(by_gene),
                "n_genes_with_multiple_records": 0,
                "pearson_r": gene_r,
                "r_squared": gene_r**2,
                "selection": "same selected set; multiple DAR records averaged per gene",
            },
        ]
    )
    regression["atac_source_path"] = str(atac_path)
    regression["rna_source_path"] = str(rna_path)
    write_tsv(regression, "figure8_atac_rna_regression_check_015.tsv")

    atac_genes = set(atac.loc[atac["log2FC"].abs() > 0.3, "geneId"])
    rna_genes = set(rna.loc[rna["logFC"].abs() > 1, "geneId"])
    overlap = len(atac_genes & rna_genes)
    universes = [22043, 23460, 24212, 55986]
    rows = []
    for universe in universes:
        expected = len(atac_genes) * len(rna_genes) / universe
        rows.append(
            {
                "universe_n": universe,
                "universe_interpretation": {
                    22043: "task009 ±1-kb promoter OCR and common-RNA universe",
                    23460: "task009 ±3-kb promoter OCR and common-RNA universe",
                    24212: "task009 jointly detectable RNA genes",
                    55986: "all raw-count annotation rows; not an appropriate joint assay universe",
                }[universe],
                "n_rna": len(rna_genes),
                "n_atac": len(atac_genes),
                "observed_overlap": overlap,
                "expected_overlap": expected,
                "fold_enrichment": overlap / expected,
                "hypergeometric_one_sided_p": hypergeom.sf(
                    overlap - 1, universe, len(atac_genes), len(rna_genes)
                ),
            }
        )
    write_tsv(pd.DataFrame(rows), "figure8_overlap_universe_sensitivity_015.tsv")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figure1()
    figure2()
    figure3()
    figure7()
    figure8()
    print("Targeted figure checks completed.")


if __name__ == "__main__":
    main()
