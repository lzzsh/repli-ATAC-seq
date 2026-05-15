#!/usr/bin/env python3
"""
Preprocess Repli-seq data for Arabidopsis and Zeamay.
Strategies aligned with the rice reference pipeline.

Label assignment per bin (matches rice GFF3 logic):
  1. CPM-normalize each sample
  2. Average replicates within each S-phase window (ES / MS / LS)
  3. log2((S_avg + 1) / (G1 + 1))
  4. If max(ES, MS, LS) score <= NON_REP_THRESHOLD → Non-replication (3)
     else label = argmax → ES(0) / MS(1) / LS(2)

Train/val/test split is handled at runtime via manifest.yaml chromosome lists.
This script outputs one TSV and one GFF3 per species+resolution (whole genome).

Label map  →  ES=0  MS=1  LS=2  Non-replication=3

Outputs per species+resolution:
  - <species>_<res>.tsv    (tab-separated, with label + score columns)
  - <species>_<res>.gff3   (GFF3 with Name= and color=, same format as
                             rice_wt_rt_1024bp.gff3)
"""

import gc
import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
DATA_DIR = Path("/Users/lzz/Downloads/qq")
OUT_DIR  = Path("/Users/lzz/Documents/GitHub/repli-ATAC-seq/Transformer/data/labels")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Column definitions ─────────────────────────────────────────────────────
COORD_COLS = ["chr", "start", "end"]

ARA_SAMPLE_COLS = [
    "G1",
    "ES_rep1", "ES_rep2", "ES_rep3",
    "MS_rep1", "MS_rep2", "MS_rep3",
    "LS_rep1", "LS_rep2", "LS_rep3",
]
ZEA_SAMPLE_COLS = [
    "G1",
    "ES_rep1", "ES_rep2", "ES_rep3",
    "MS_rep1", "MS_rep2", "MS_rep3",
    "LS_rep1", "LS_rep2",            # zeamay has no LS_rep3
]

PHASE_REPS = {
    "ES": ["ES_rep1", "ES_rep2", "ES_rep3"],
    "MS": ["MS_rep1", "MS_rep2", "MS_rep3"],
    "LS": ["LS_rep1", "LS_rep2", "LS_rep3"],
}

# Bins whose best S-phase log2-ratio is below this are Non-replication.
# Mirrors rice: bins with negligible S-phase enrichment over G1.
NON_REP_THRESHOLD = 0.0   # log2(S/G1) ≤ 0 → not enriched in any S window

LABEL_NAMES  = {0: "ES", 1: "MS", 2: "LS", 3: "Non-replication"}
LABEL_COLORS = {
    "ES":              "#2250F1",
    "MS":              "#9B30FF",
    "LS":              "#FB0018",
    "Non-replication": "#B0B0B0",
}

# ── Chromosomes to keep (exclude organelles / scaffolds) ──────────────────
ARA_CHROMS = ["1", "2", "3", "4", "5"]        # exclude Mt, Pt
ZEA_CHROMS = [str(i) for i in range(1, 11)]   # exclude scaffolds


# ── Core helpers ───────────────────────────────────────────────────────────

def cpm(arr: np.ndarray) -> np.ndarray:
    lib = arr.sum(axis=0, keepdims=True).clip(min=1)
    return arr / lib * 1e6


def label_chunk(df: pd.DataFrame, sample_cols: list,
                cpm_counts: np.ndarray | None = None) -> None:
    """Add label + score columns to df in-place (no rows dropped).

    cpm_counts: pre-computed whole-genome CPM array aligned to df's rows.
    If None, CPM is computed from df alone (per-chunk, legacy behaviour).
    """
    if cpm_counts is not None:
        counts = cpm_counts
    else:
        counts = cpm(df[sample_cols].values.astype(np.float32))
    col_idx = {c: i for i, c in enumerate(sample_cols)}
    g1      = counts[:, col_idx["G1"]]

    scores = []
    for phase in ("ES", "MS", "LS"):
        present = [r for r in PHASE_REPS[phase] if r in col_idx]
        idx     = [col_idx[r] for r in present]
        avg     = counts[:, idx].mean(axis=1)
        scores.append(np.log2((avg + 1.0) / (g1 + 1.0)))

    score_mat = np.stack(scores, axis=1)          # (N, 3)
    best      = score_mat.max(axis=1)

    labels = np.where(
        best <= NON_REP_THRESHOLD,
        np.int8(3),                               # Non-replication
        score_mat.argmax(axis=1).astype(np.int8), # ES / MS / LS
    )

    df["label"]    = labels
    df["ES_score"] = score_mat[:, 0]
    df["MS_score"] = score_mat[:, 1]
    df["LS_score"] = score_mat[:, 2]


# ── GFF3 writer ───────────────────────────────────────────────────────────

def write_gff3(df: pd.DataFrame, out_path: Path) -> None:
    """
    Write a GFF3 matching the rice format:
      chr  assign_rt_class  peaks  start1  end1  .  .  .  Name=X;color=#XXXXXX;
    Coordinates are 1-based closed intervals (BED start+1 .. end).
    """
    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for row in df.itertuples(index=False):
            name  = LABEL_NAMES[row.label]
            color = LABEL_COLORS[name]
            # BED is 0-based half-open → GFF3 is 1-based closed
            fh.write(
                f"{row.chr}\tassign_rt_class\tpeaks\t"
                f"{row.start + 1}\t{row.end}\t.\t.\t.\t"
                f"Name={name};color={color};\n"
            )




def process(species: str, bed_128: Path, bed_1024: Path,
            sample_cols: list, keep_chroms: list):

    print(f"\n{'='*60}\n  {species}\n{'='*60}")

    for res, bed_path in [("128bp", bed_128), ("1024bp", bed_1024)]:
        print(f"\n  [{res}]  {bed_path.name}")
        col_names = COORD_COLS + sample_cols

        print("    reading file …", flush=True)
        df_all = pd.read_csv(
            bed_path, sep="\t", header=None, names=col_names,
            dtype={"chr": str}, low_memory=False,
        )
        df_all = df_all[df_all["chr"].isin(keep_chroms)].reset_index(drop=True)
        print(f"    {len(df_all):,} bins after chrom filter", flush=True)

        # whole-genome CPM normalisation (library size = sum over all bins)
        wg_cpm = cpm(df_all[sample_cols].values.astype(np.float32))

        frames = []
        for chrom, chunk in df_all.groupby("chr", sort=False):
            idx = chunk.index
            chunk = chunk.reset_index(drop=True)
            label_chunk(chunk, sample_cols, cpm_counts=wg_cpm[idx])
            frames.append(chunk)

        del df_all
        gc.collect()

        out_df   = pd.concat(frames, ignore_index=True)
        tsv_path = OUT_DIR / f"{species}_{res}.tsv"
        out_df.to_csv(tsv_path, sep="\t", index=False)

        n  = len(out_df)
        es = (out_df["label"] == 0).sum()
        ms = (out_df["label"] == 1).sum()
        ls = (out_df["label"] == 2).sum()
        nr = (out_df["label"] == 3).sum()
        print(f"    {tsv_path.name}: {n:>9,} bins  "
              f"ES={es:,}  MS={ms:,}  LS={ls:,}  Non-rep={nr:,}")

        gff3_path = OUT_DIR / f"{species}_{res}.gff3"
        write_gff3(out_df, gff3_path)
        print(f"    {gff3_path.name}: {n:,} bins written")

        del out_df
        gc.collect()


# ── Run ────────────────────────────────────────────────────────────────────

process(
    species     = "arabidopsis",
    bed_128     = DATA_DIR / "arabidopsis_128bp_matrix_sorted.bed",
    bed_1024    = DATA_DIR / "arabidopsis_1024bp_matrix_sorted.bed",
    sample_cols = ARA_SAMPLE_COLS,
    keep_chroms = ARA_CHROMS,
)

process(
    species     = "zeamay",
    bed_128     = DATA_DIR / "zeamay2017_128bp_matrix_sorted.bed",
    bed_1024    = DATA_DIR / "zeamay2017_1024bp_matrix_sorted.bed",
    sample_cols = ZEA_SAMPLE_COLS,
    keep_chroms = ZEA_CHROMS,
)

print("\nAll done. Output →", OUT_DIR)
