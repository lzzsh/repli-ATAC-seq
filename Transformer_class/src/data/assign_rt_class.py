"""
Assign replication timing (RT) class to each bin from raw count TSV,
and write results as GFF3.

Logic:
  1. TPM-normalize ES, MS, LS, G1 counts independently across all bins.
  2. For each bin, compute ratio = phase_tpm / g1_tpm for ES, MS, LS.
  3. If all three ratios < 1.0 → NR (not replicated).
  4. Otherwise, the phase with the highest ratio is the RT class (ES / MS / LS).

Input TSV columns: chrom, start, end, ES_count, MS_count, LS_count, G1_count
Output: GFF3 file (one feature per bin, Name= and color= attributes for IGV)
"""

import argparse
import numpy as np
import pandas as pd


_COLOR = {
    "ES": "#2250F1",
    "MS": "#1A8A12",
    "LS": "#FB0018",
    "NR": "#B0B0B0",
}

_IGV_NAME = {"NR": "Non-replication"}


def tpm(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts, dtype=np.float64)
    return counts / total * 1e6


def assign_rt_class(df: pd.DataFrame) -> pd.Series:
    es = tpm(df["ES_count"].values.astype(float))
    ms = tpm(df["MS_count"].values.astype(float))
    ls = tpm(df["LS_count"].values.astype(float))
    g1 = tpm(df["G1_count"].values.astype(float))

    eps = 1e-6
    ratios = np.stack([es, ms, ls], axis=1) / (g1[:, None] + eps)  # [N, 3]
    max_ratio = ratios.max(axis=1)
    best = ratios.argmax(axis=1)

    phase_names = ["ES", "MS", "LS"]
    classes = [
        "NR" if max_ratio[i] < 1.0 else phase_names[best[i]]
        for i in range(len(df))
    ]
    return pd.Series(classes, index=df.index, name="RT_class")


def write_gff3(df: pd.DataFrame, path: str, source: str = "assign_rt_class") -> None:
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for row in df.itertuples(index=False):
            name = _IGV_NAME.get(row.RT_class, row.RT_class)
            color = _COLOR[row.RT_class]
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
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if "#chrom" in df.columns:
        df = df.rename(columns={"#chrom": "chrom"})

    for col in ["chrom", "start", "end", "ES_count", "MS_count", "LS_count", "G1_count"]:
        assert col in df.columns, f"Missing column: {col}"

    df["RT_class"] = assign_rt_class(df)

    counts = df["RT_class"].value_counts()
    print("RT class distribution:")
    for cls in ["ES", "MS", "LS", "NR"]:
        n = counts.get(cls, 0)
        print(f"  {cls}: {n:>8,}  ({n/len(df)*100:.1f}%)")

    write_gff3(df, args.output)
    print(f"\nWritten {len(df):,} features to {args.output}")


if __name__ == "__main__":
    main()
