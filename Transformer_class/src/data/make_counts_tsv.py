"""
Convert repli_peaks_quan.txt to the count TSV expected by load_labels.

Input columns (1-based):
  1: chr  2: start  3: end
  4: G1-1  5: G1-2
  6: ES-1  7: ES-2
  8: MS
  9: LS

Output TSV: chrom  start  end  ES_count  MS_count  LS_count  G1_count
  - G1_count = mean(G1-1, G1-2)
  - ES_count = mean(ES-1, ES-2)
  - MS_count and LS_count are single replicates
"""

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="repli_peaks_quan.txt")
    parser.add_argument("output", help="output TSV path")
    args = parser.parse_args()

    df = pd.read_csv(
        args.input, sep="\t", header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8],
        names=["chrom", "start", "end", "G1_1", "G1_2", "ES_1", "ES_2", "MS", "LS"],
    )

    out = pd.DataFrame({
        "chrom":    df["chrom"],
        "start":    df["start"],
        "end":      df["end"],
        "ES_count": (df["ES_1"] + df["ES_2"]) / 2.0,
        "MS_count": df["MS"].astype(float),
        "LS_count": df["LS"].astype(float),
        "G1_count": (df["G1_1"] + df["G1_2"]) / 2.0,
    })

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(out):,} bins to {args.output}")


if __name__ == "__main__":
    main()
