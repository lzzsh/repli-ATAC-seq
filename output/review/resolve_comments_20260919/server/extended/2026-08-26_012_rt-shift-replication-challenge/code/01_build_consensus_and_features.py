#!/usr/bin/env python3
"""Build conservative CPP8 consensus RT instability and genomic features per 1-kb bin."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


HISTONE_FILES = {
    "h3k27ac_fraction": "H3K27ac_roots.bed",
    "h3k4me1_fraction": "H3K4me1_roots.bed",
    "h3k4me3_fraction": "H3K4me3_roots.bed",
    "h3k9me2_fraction": "H3K9me2_roots.bed",
    "h3k27me3_fraction": "H3K27me3_roots.bed",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--rt-bins", required=True)
    p.add_argument("--raw-repli", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--annotation-bed", required=True)
    p.add_argument("--histone-dir", required=True)
    p.add_argument("--ocr-metrics", required=True)
    p.add_argument("--rna", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def normalize_chrom(chrom: str) -> str:
    m = re.fullmatch(r"chr0?([1-9]|1[0-2])", str(chrom))
    return f"chr{int(m.group(1)):02d}" if m else str(chrom)


def parse_genes(path: str) -> pd.DataFrame:
    names = ["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    gff = pd.read_csv(path, sep="\t", comment="#", header=None, names=names, low_memory=False)
    genes = gff[(gff["feature"] == "gene") & gff["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy()
    genes["gene_id"] = genes["attributes"].astype(str).str.extract(r"(?:^|;)ID=([^;]+)", expand=False)
    genes[["start", "end"]] = genes[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    genes = genes.dropna(subset=["gene_id", "start", "end"])
    genes[["start", "end"]] = genes[["start", "end"]].astype(int)
    genes["start0"] = genes["start"] - 1
    genes["length"] = genes["end"] - genes["start0"]
    genes["midpoint"] = (genes["start0"] + genes["end"]) // 2
    genes["is_centromere_specific"] = genes["attributes"].astype(str).str.contains("centromere-specific", case=False, regex=False)
    genes = genes.drop_duplicates("gene_id")
    return genes[["gene_id", "chr", "start0", "end", "length", "midpoint", "is_centromere_specific"]]


def interval_coverage(bins: pd.DataFrame, intervals: pd.DataFrame) -> np.ndarray:
    coverage = np.zeros(len(bins), dtype=float)
    for chrom, bi in bins.groupby("chr", sort=False).groups.items():
        bidx = np.asarray(bi, dtype=int)
        bd = bins.loc[bidx]
        starts = bd["start"].to_numpy(int)
        ends = bd["end"].to_numpy(int)
        iv = intervals[intervals["chr"] == chrom]
        for row in iv[["start", "end"]].itertuples(index=False):
            left = int(np.searchsorted(ends, int(row.start), side="right"))
            right = int(np.searchsorted(starts, int(row.end), side="left"))
            if right <= left:
                continue
            ov = np.minimum(ends[left:right], int(row.end)) - np.maximum(starts[left:right], int(row.start))
            coverage[bidx[left:right]] += np.maximum(ov, 0)
    length = (bins["end"] - bins["start"]).to_numpy(float)
    return np.clip(coverage / length, 0, 1)


def local_count(query: np.ndarray, positions: np.ndarray, radius: int) -> np.ndarray:
    left = np.searchsorted(positions, query - radius, side="left")
    right = np.searchsorted(positions, query + radius, side="right")
    return right - left


def local_mean(query: np.ndarray, positions: np.ndarray, values: np.ndarray, radius: int) -> np.ndarray:
    valid = np.isfinite(values)
    sums = np.concatenate([[0.0], np.cumsum(np.where(valid, values, 0.0))])
    counts = np.concatenate([[0], np.cumsum(valid.astype(int))])
    left = np.searchsorted(positions, query - radius, side="left")
    right = np.searchsorted(positions, query + radius, side="right")
    total = sums[right] - sums[left]
    n = counts[right] - counts[left]
    return np.divide(total, n, out=np.full(len(query), np.nan), where=n > 0)


def nearest_values(query: np.ndarray, positions: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    idx = np.searchsorted(positions, query)
    li = np.maximum(idx - 1, 0)
    ri = np.minimum(idx, len(positions) - 1)
    choose_right = np.abs(query - positions[ri]) < np.abs(query - positions[li])
    chosen = np.where(choose_right, ri, li)
    return values[chosen], np.abs(query - positions[chosen])


def read_bed(path: str, label: str | None = None) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3] if label is not None else [0, 1, 2])
    bed = bed.rename(columns={0: "chr", 1: "start", 2: "end", 3: "label"})
    bed["chr"] = bed["chr"].astype(str).map(normalize_chrom)
    if label is not None:
        bed = bed[bed["label"] == label]
    bed[["start", "end"]] = bed[["start", "end"]].astype(int)
    return bed[["chr", "start", "end"]]


def raw_wt_g1(path: str) -> pd.DataFrame:
    names = [
        "chr", "start", "end", "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
        "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
        "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
        "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
        "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS",
    ]
    d = pd.read_csv(path, sep="\t", header=None, names=names)
    d["wt_g1_mean_count"] = d[["WT-1-G1", "WT-2-G1"]].mean(axis=1)
    return d[["chr", "start", "end", "wt_g1_mean_count"]]


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in outdir.parents and outdir != task_root:
        raise ValueError(f"Refusing output outside task root: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    rt = pd.read_csv(args.rt_bins, sep="\t")
    rt = rt[rt["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy().reset_index(drop=True)
    d1 = pd.to_numeric(rt["delta_rt_cpp8_1"], errors="coerce").to_numpy(float)
    d3 = pd.to_numeric(rt["delta_rt_cpp8_3"], errors="coerce").to_numpy(float)
    same = np.sign(d1) == np.sign(d3)
    consensus = np.where(same, np.sign((d1 + d3) / 2) * np.minimum(np.abs(d1), np.abs(d3)), 0.0)
    out = rt[["chr", "start", "end", "block", "WT_wrt", "delta_rt_cpp8_1", "delta_rt_cpp8_3"]].copy()
    out["midpoint"] = ((out["start"] + out["end"]) // 2).astype(int)
    out["consensus_shift"] = consensus
    out["consensus_instability"] = np.abs(consensus)
    out["alleles_same_sign"] = same
    for pct in (2.5, 5, 10):
        threshold = float(np.nanpercentile(out["consensus_instability"], 100 - pct))
        out[f"hotspot_top_{str(pct).replace('.', '_')}pct"] = out["consensus_instability"] >= threshold
    primary_threshold = float(np.nanpercentile(out["consensus_instability"], 95))
    out["signed_shift_class"] = np.where(~same, "discordant", np.where(out["consensus_instability"] < primary_threshold, "stable", np.where(consensus < 0, "concordant_earlier", "concordant_later")))

    g1 = raw_wt_g1(args.raw_repli)
    out = out.merge(g1, on=["chr", "start", "end"], how="left")
    out["log_wt_g1_coverage"] = np.log1p(out["wt_g1_mean_count"])

    fasta = pysam.FastaFile(args.fasta)
    gc = np.full(len(out), np.nan)
    for chrom, bi in out.groupby("chr", sort=False).groups.items():
        idx = np.asarray(bi, dtype=int)
        seq = np.frombuffer(fasta.fetch(chrom).upper().encode("ascii"), dtype=np.uint8)
        valid = np.isin(seq, np.frombuffer(b"ACGT", dtype=np.uint8)).astype(int)
        is_gc = np.isin(seq, np.frombuffer(b"GC", dtype=np.uint8)).astype(int)
        valid_prefix = np.concatenate([[0], np.cumsum(valid)])
        gc_prefix = np.concatenate([[0], np.cumsum(is_gc)])
        starts = out.loc[idx, "start"].to_numpy(int)
        ends = out.loc[idx, "end"].to_numpy(int)
        n_valid = valid_prefix[ends] - valid_prefix[starts]
        n_gc = gc_prefix[ends] - gc_prefix[starts]
        gc[idx] = np.divide(n_gc, n_valid, out=np.full(len(idx), np.nan), where=n_valid > 0)
    fasta.close()
    out["gc_fraction"] = gc

    genes = parse_genes(args.gff)
    gene_intervals = genes.rename(columns={"start0": "start"})[["chr", "start", "end"]]
    out["gene_body_fraction"] = interval_coverage(out, gene_intervals)
    rna = pd.read_csv(args.rna, sep="\t")[["gene_id", "wt_mean_cpm"]]
    genes = genes.merge(rna, on="gene_id", how="left")
    ocr = pd.read_csv(args.ocr_metrics, sep="\t", usecols=["chr", "start", "end"])
    ocr["midpoint"] = ((ocr["start"] + ocr["end"]) // 2).astype(int)
    for chrom, bi in out.groupby("chr", sort=False).groups.items():
        idx = np.asarray(bi, dtype=int)
        query = out.loc[idx, "midpoint"].to_numpy(int)
        gd = genes[genes["chr"] == chrom].sort_values("midpoint")
        gp = gd["midpoint"].to_numpy(int)
        out.loc[idx, "gene_count_50kb"] = local_count(query, gp, 50_000)
        nearest_length, gene_dist = nearest_values(query, gp, gd["length"].to_numpy(float))
        out.loc[idx, "nearest_gene_length"] = nearest_length
        out.loc[idx, "distance_to_nearest_gene"] = gene_dist
        out.loc[idx, "wt_expression_50kb_log2cpm"] = local_mean(query, gp, np.log2(gd["wt_mean_cpm"].to_numpy(float) + 1.0), 50_000)
        op = np.sort(ocr.loc[ocr["chr"] == chrom, "midpoint"].to_numpy(int))
        out.loc[idx, "ocr_count_50kb"] = local_count(query, op, 50_000)
        cent = gd.loc[gd["is_centromere_specific"], "midpoint"].to_numpy(int)
        if len(cent):
            proxy = int(np.median(cent))
            out.loc[idx, "centromere_repeat_cluster_proxy"] = proxy
            out.loc[idx, "distance_to_centromere_proxy"] = np.abs(query - proxy)

    annotation = pd.read_csv(args.annotation_bed, sep="\t", header=None, names=["chr", "start", "end", "label"])
    annotation["chr"] = annotation["chr"].astype(str).map(normalize_chrom)
    for label, col in [("Class I elements", "class_i_te_fraction"), ("Class II elements", "class_ii_te_fraction")]:
        out[col] = interval_coverage(out, annotation[annotation["label"] == label][["chr", "start", "end"]])

    histone_dir = Path(args.histone_dir)
    for col, filename in HISTONE_FILES.items():
        out[col] = interval_coverage(out, read_bed(str(histone_dir / filename)))

    out["active_histone_fraction"] = out[["h3k27ac_fraction", "h3k4me1_fraction", "h3k4me3_fraction"]].mean(axis=1)
    out["repressive_histone_fraction"] = out[["h3k9me2_fraction", "h3k27me3_fraction"]].mean(axis=1)
    out.to_csv(outdir / "cpp8_consensus_rt_instability_012.tsv.gz", sep="\t", index=False)
    out.to_csv(outdir / "replication_challenge_features_012.tsv.gz", sep="\t", index=False)
    genes.to_csv(outdir / "gene_features_source_012.tsv.gz", sep="\t", index=False)
    print("Bins:", len(out))
    print("Primary top-5% threshold:", primary_threshold)
    print(out["signed_shift_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
