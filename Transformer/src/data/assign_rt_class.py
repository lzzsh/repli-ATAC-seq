"""
Assign replication timing (RT) class to each bin from raw count TSV,
and write results as GFF3.

Logic:
  1. TPM-normalize ES, MS, LS, G1 counts independently across all bins.
  2. For each bin, compute ratio = phase_tpm / g1_tpm for ES, MS, LS.
  3. If all three ratios < 1.0  → NR (not replicated).
  4. Otherwise, find max ratio among ES/MS/LS.
     Phases whose ratio >= 0.9 * max_ratio are "dominant".
  5. Single dominant phase → that class (ES / MS / LS).
     Multiple dominant phases → combined class (ESMS / MSLS / ESLS / ESMSLS).
     These combined classes are treated as boundary bins (filtered during training).

Input TSV columns: chrom, start, end, ES_count, MS_count, LS_count, G1_count
Output: GFF3 file (one feature per bin, Name= attribute holds the RT class)
"""

import argparse
import numpy as np
import pandas as pd


_COMBO = {
    (True, False, False): "ES",
    (False, True, False): "MS",
    (False, False, True): "LS",
    (True, True, False):  "ESMS",
    (False, True, True):  "MSLS",
    (True, False, True):  "ESLS",
    (True, True, True):   "ESMSLS",
}


def tpm(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts, dtype=np.float64)
    return counts / total * 1e6


def assign_rt_class(df: pd.DataFrame, threshold: float = 0.9) -> pd.Series:
    es = tpm(df["ES_count"].values.astype(float))
    ms = tpm(df["MS_count"].values.astype(float))
    ls = tpm(df["LS_count"].values.astype(float))
    g1 = tpm(df["G1_count"].values.astype(float))

    eps = 1e-6
    ratio_es = es / (g1 + eps)
    ratio_ms = ms / (g1 + eps)
    ratio_ls = ls / (g1 + eps)

    classes = []
    for re, rm, rl in zip(ratio_es, ratio_ms, ratio_ls):
        max_r = max(re, rm, rl)
        if max_r < 1.0:
            classes.append("NR")
            continue
        dominant = (
            re >= threshold * max_r,
            rm >= threshold * max_r,
            rl >= threshold * max_r,
        )
        classes.append(_COMBO.get(dominant, "NR"))

    return pd.Series(classes, index=df.index, name="RT_class")


_COLOR = {
    "ES":      "#2250F1",
    "MS":      "#1A8A12",
    "LS":      "#FB0018",
    "NR":      "#B0B0B0",
    "ESMS":    "#28C5CC",
    "MSLS":    "#FFFD33",
    "ESLS":    "#EA3CF2",
    "ESMSLS":  "#FAB427",
}

_IGV_NAME = {
    "NR": "Non-replication",
}


def write_gff3(df: pd.DataFrame, path: str, source: str = "assign_rt_class") -> None:
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, row in enumerate(df.itertuples(index=False), start=1):
            name = _IGV_NAME.get(row.RT_class, row.RT_class)
            color = _COLOR.get(row.RT_class, "#B0B0B0")
            # GFF3 coords are 1-based closed; input TSV is 0-based half-open
            f.write(
                f"{row.chrom}\t{source}\tpeaks\t"
                f"{row.start + 1}\t{row.end}\t"
                f".\t.\t.\t"
                f"Name={name};color={color};\n"
            )


def main():
    parser = argparse.ArgumentParser(description="Assign RT class from raw count TSV, output GFF3.")
    parser.add_argument("input",  help="Input TSV with ES/MS/LS/G1 counts")
    parser.add_argument("output", help="Output GFF3 path")
    parser.add_argument("--threshold", type=float, default=0.9,
                        help="Fraction of max ratio to consider dominant (default: 0.9)")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if "#chrom" in df.columns:
        df = df.rename(columns={"#chrom": "chrom"})

    for col in ["chrom", "start", "end", "ES_count", "MS_count", "LS_count", "G1_count"]:
        assert col in df.columns, f"Missing column: {col}"

    df["RT_class"] = assign_rt_class(df, threshold=args.threshold)

    counts = df["RT_class"].value_counts()
    print("RT class distribution:")
    for cls in ["ES", "MS", "LS", "ESMS", "MSLS", "ESLS", "ESMSLS", "NR"]:
        n = counts.get(cls, 0)
        print(f"  {cls:8s}: {n:>8,}  ({n/len(df)*100:.1f}%)")

    write_gff3(df, args.output)
    print(f"\nWritten {len(df):,} features to {args.output}")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
