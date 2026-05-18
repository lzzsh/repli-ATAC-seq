"""
Convert raw-count BED files (128 bp bins) to model-ready TSV signals.

Pipeline per species (Enformer style):
  1. CPM normalise each sample column (counts per million reads)
  2. Soft clip each replicate: x > CLIP → CLIP - 1 + sqrt(x - CLIP + 1)
  3. Average replicates → G1, ES, MS, LS
  4. log1p transform: log(1 + x)
  5. Write TSV: chrom  start  end  G1  ES  MS  LS

Output values are non-negative (log1p of clipped CPM).
Model output uses Softplus activation to match this target space.
Each species is normalised independently; cross-species comparability
is handled by per-species output heads, not by the target values.

Usage:
  python preprocess_signals.py          # processes all three species
  python preprocess_signals.py rice     # single species
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

OUT_DIR   = Path(__file__).parent / "labels"
CLIP_SOFT = 384.0   # soft clip threshold on CPM values (matches basenji2)

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


def cpm_normalise(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    for col in sample_cols:
        total = df[col].sum()
        if total > 0:
            df[col] = df[col] / total * 1e6
    return df


def soft_clip(arr: np.ndarray, threshold: float) -> np.ndarray:
    out = arr.copy()
    mask = out > threshold
    out[mask] = threshold - 1.0 + np.sqrt(arr[mask] - threshold + 1.0)
    return out


def process(name: str, cfg: dict) -> None:
    print(f"[{name}] reading {cfg['input'].name} ...", flush=True)
    df = pd.read_csv(cfg["input"], sep="\t", header=None,
                     names=["chrom", "start", "end"] + cfg["cols"],
                     dtype={"chrom": str}, low_memory=False)

    print(f"[{name}] CPM normalising {len(cfg['cols'])} samples ...", flush=True)
    df = cpm_normalise(df, cfg["cols"])

    print(f"[{name}] soft clip each replicate (threshold={CLIP_SOFT}) ...", flush=True)
    for col in cfg["cols"]:
        df[col] = soft_clip(df[col].values.astype(np.float64), CLIP_SOFT)

    print(f"[{name}] averaging replicates ...", flush=True)
    out = df[["chrom", "start", "end"]].copy()
    for phase, reps in cfg["groups"].items():
        out[phase] = df[reps].mean(axis=1).values

    print(f"[{name}] log1p ...", flush=True)
    for phase in cfg["groups"]:
        vals = np.log1p(out[phase].values.astype(np.float64))
        out[phase] = vals.astype(np.float32)
        print(f"  {phase}: mean={out[phase].mean():.4f}  std={out[phase].std():.4f}"
              f"  min={out[phase].min():.4f}  max={out[phase].max():.4f}")

    out.to_csv(cfg["output"], sep="\t", index=False, float_format="%.6f")
    print(f"[{name}] wrote {len(out):,} rows → {cfg['output'].name}\n")


def main() -> None:
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(SPECIES)
    for name in targets:
        if name not in SPECIES:
            print(f"unknown species '{name}', choices: {list(SPECIES)}")
            sys.exit(1)
        process(name, SPECIES[name])


if __name__ == "__main__":
    main()
