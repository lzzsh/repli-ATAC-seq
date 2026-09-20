#!/usr/bin/env python3
"""Build nested-window, annulus and offset RT metrics around each OCR."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ALLELES = ("cpp8_1", "cpp8_3")
SCALES_KB = (1, 5, 25, 100)
RADII_KB = (1, 5, 20, 50, 100)
OFFSETS_KB = tuple(range(-100, 101, 5))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ocr-metrics", required=True)
    p.add_argument("--rt-bins", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def prefix_for(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(values)
    sums = np.concatenate([[0.0], np.cumsum(np.where(valid, values, 0.0))])
    counts = np.concatenate([[0], np.cumsum(valid.astype(int))])
    return sums, counts


def range_mean(centers: np.ndarray, sums: np.ndarray, counts: np.ndarray, query: np.ndarray, left_delta: int, right_delta: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    left = np.searchsorted(centers, query + left_delta, side="left")
    right = np.searchsorted(centers, query + right_delta, side="right")
    total = sums[right] - sums[left]
    n = counts[right] - counts[left]
    mean = np.divide(total, n, out=np.full(len(query), np.nan), where=n > 0)
    return mean, total, n


def safe_annulus(outer_sum: np.ndarray, outer_n: np.ndarray, inner_sum: np.ndarray, inner_n: np.ndarray) -> np.ndarray:
    total = outer_sum - inner_sum
    n = outer_n - inner_n
    return np.divide(total, n, out=np.full(len(n), np.nan), where=n > 0)


def main() -> None:
    args = parse_args()
    output = Path(args.output).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in output.parents:
        raise ValueError(f"Refusing output outside task root: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    ocr = pd.read_csv(args.ocr_metrics, sep="\t")
    rt = pd.read_csv(args.rt_bins, sep="\t")
    keep_cols = [
        "chr", "start", "end", "block", "length", "WT_wrt",
        "delta_rt_cpp8_1", "delta_rt_cpp8_3",
        "wt_es_log_atac", "wt_ms_log_atac", "wt_ls_log_atac",
        "delta_atac_cpp8_1_es", "delta_atac_cpp8_1_ms", "delta_atac_cpp8_1_ls",
        "delta_atac_cpp8_3_es", "delta_atac_cpp8_3_ms", "delta_atac_cpp8_3_ls",
    ]
    missing = set(keep_cols) - set(ocr.columns)
    if missing:
        raise ValueError(f"Missing OCR columns: {sorted(missing)}")
    out = ocr[keep_cols].copy()
    out["midpoint"] = ((out["start"] + out["end"]) // 2).astype(int)
    value_cols = ["WT_wrt", "delta_rt_cpp8_1", "delta_rt_cpp8_3"]

    feature_data: dict[str, np.ndarray] = {}
    for chrom, oi in out.groupby("chr", sort=False).groups.items():
        qi = np.asarray(oi, dtype=int)
        query = out.loc[qi, "midpoint"].to_numpy(int)
        rd = rt[rt["chr"] == chrom].sort_values("start")
        centers = ((rd["start"].to_numpy(int) + rd["end"].to_numpy(int)) // 2).astype(int)
        for value_col in value_cols:
            values = pd.to_numeric(rd[value_col], errors="coerce").to_numpy(float)
            sums, counts = prefix_for(values)
            nested_cache: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
            for radius_kb in sorted(set(SCALES_KB) | set(RADII_KB)):
                radius = radius_kb * 1000
                nested_cache[radius_kb] = range_mean(centers, sums, counts, query, -radius, radius)
            for scale_kb in SCALES_KB:
                name = f"{value_col}_scale_{scale_kb}kb"
                feature_data.setdefault(name, np.full(len(out), np.nan))[qi] = nested_cache[scale_kb][0]
            annuli = [(1, 5), (5, 20), (20, 50), (50, 100)]
            for inner, outer in annuli:
                _, osum, on = nested_cache[outer]
                _, isum, inn = nested_cache[inner]
                name = f"{value_col}_annulus_{inner}_{outer}kb"
                feature_data.setdefault(name, np.full(len(out), np.nan))[qi] = safe_annulus(osum, on, isum, inn)
            for offset_kb in OFFSETS_KB:
                center_shift = offset_kb * 1000
                mean, _, _ = range_mean(centers, sums, counts, query, center_shift - 2500, center_shift + 2500)
                tag = f"m{abs(offset_kb)}" if offset_kb < 0 else f"p{offset_kb}"
                name = f"{value_col}_offset_{tag}kb"
                feature_data.setdefault(name, np.full(len(out), np.nan))[qi] = mean

    for name, values in feature_data.items():
        out[name] = values
    out.to_csv(output, sep="\t", index=False)
    print("OCRs:", len(out))
    print("Columns:", len(out.columns))
    print("Nested scales (radius, kb):", SCALES_KB)
    print("Offsets (5-kb window centres, kb):", OFFSETS_KB)


if __name__ == "__main__":
    main()
