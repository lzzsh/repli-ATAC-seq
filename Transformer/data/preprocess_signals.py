"""
Convert raw-count BED files (128 bp bins) to model-ready TSV signals.

Pipeline per species:
  1. TPM normalise each sample column (bin_size = 128 bp)
  2. Average replicates → G1, ES, MS, LS
  3. Write TSV: chrom  start  end  G1  ES  MS  LS

Usage:
  python preprocess_signals.py          # processes all three species
  python preprocess_signals.py rice     # single species
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

BIN_SIZE = 128
OUT_DIR = Path(__file__).parent / "labels"

# ── species configs ───────────────────────────────────────────────────────────
SPECIES = {
    "rice": {
        "input": OUT_DIR / "repli_peaks_quan_128.txt",
        "output": OUT_DIR / "rice_128bp_rt_signals.tsv",
        "cols": ["G1_rep1", "G1_rep2", "ES_rep1", "ES_rep2", "MS", "LS"],
        "groups": {
            "G1": ["G1_rep1", "G1_rep2"],
            "ES": ["ES_rep1", "ES_rep2"],
            "MS": ["MS"],
            "LS": ["LS"],
        },
    },
    "arabidopsis": {
        "input": OUT_DIR / "arabidopsis_128bp_matrix_sorted.bed",
        "output": OUT_DIR / "arabidopsis_128bp_rt_signals.tsv",
        "cols": ["G1", "ES_rep1", "ES_rep2", "ES_rep3",
                 "MS_rep1", "MS_rep2", "MS_rep3",
                 "LS_rep1", "LS_rep2", "LS_rep3"],
        "groups": {
            "G1": ["G1"],
            "ES": ["ES_rep1", "ES_rep2", "ES_rep3"],
            "MS": ["MS_rep1", "MS_rep2", "MS_rep3"],
            "LS": ["LS_rep1", "LS_rep2", "LS_rep3"],
        },
    },
    "zeamay": {
        "input": OUT_DIR / "zeamay2017_128bp_matrix_sorted.bed",
        "output": OUT_DIR / "zeamay_128bp_rt_signals.tsv",
        "cols": ["G1", "ES_rep1", "ES_rep2", "ES_rep3",
                 "MS_rep1", "MS_rep2", "MS_rep3",
                 "LS_rep1", "LS_rep2"],
        "groups": {
            "G1": ["G1"],
            "ES": ["ES_rep1", "ES_rep2", "ES_rep3"],
            "MS": ["MS_rep1", "MS_rep2", "MS_rep3"],
            "LS": ["LS_rep1", "LS_rep2"],
        },
    },
}


def tpm_normalise(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    rpk_per_kb = 1000.0 / BIN_SIZE  # constant for fixed-size bins
    for col in sample_cols:
        rpk = df[col] * rpk_per_kb
        scale = rpk.sum() / 1e6
        df[col] = rpk / scale if scale > 0 else rpk
    return df


def process(name: str, cfg: dict) -> None:
    print(f"[{name}] reading {cfg['input'].name} ...", flush=True)
    df = pd.read_csv(cfg["input"], sep="\t", header=None,
                     names=["chrom", "start", "end"] + cfg["cols"],
                     dtype={"chrom": str}, low_memory=False)

    print(f"[{name}] TPM normalising {len(cfg['cols'])} samples ...", flush=True)
    df = tpm_normalise(df, cfg["cols"])

    print(f"[{name}] averaging replicates ...", flush=True)
    out = df[["chrom", "start", "end"]].copy()
    for phase, reps in cfg["groups"].items():
        out[phase] = df[reps].mean(axis=1)

    out.to_csv(cfg["output"], sep="\t", index=False, float_format="%.6f")
    print(f"[{name}] wrote {len(out):,} rows → {cfg['output'].name}")


def main() -> None:
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(SPECIES)
    for name in targets:
        if name not in SPECIES:
            print(f"unknown species '{name}', choices: {list(SPECIES)}")
            sys.exit(1)
        process(name, SPECIES[name])


if __name__ == "__main__":
    main()
