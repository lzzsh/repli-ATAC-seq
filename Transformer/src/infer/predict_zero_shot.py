#!/usr/bin/env python3
"""
Zero-shot cross-species prediction using the rice-trained model.

Uses the rice head to predict on arabidopsis / zeamay sequences.
The trunk is shared and species-agnostic; only the output head differs —
here we reuse the rice head as-is (zero-shot transfer).

Usage:
  python predict_zero_shot.py \
      --checkpoint /path/to/best_model.pt \
      --config     src/configs/transformer_wt.yaml \
      --fasta      /path/to/arabidopsis.fa \
      --species    arabidopsis \
      --chroms     1 2 3 4 5 \
      --out_dir    output/arabidopsis_zeroshot \
      [--batch_size 8] [--no_rc]
"""

import argparse
import sys
import numpy as np
import torch
import pyBigWig
from pathlib import Path
from pyfaidx import Fasta

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.infer._utils import load_model, fetch_one_hot, rc_average
from src.data.data_utils import GenomeSequence

BIN_SIZE  = 128
CROP_BINS = 320
OUT_BINS  = 896
STRIDE    = OUT_BINS * BIN_SIZE   # 114688 bp, gapless


def predict_genome(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print(f"Loading model from {args.checkpoint} ...")
    model, cfg = load_model(args.checkpoint, args.config, device)
    window_size = cfg["data"]["input_window_length"]
    print(f"Window size: {window_size} bp  |  Using head: rice (zero-shot)")

    genome = GenomeSequence(args.fasta)
    fa = Fasta(args.fasta)
    all_chroms = list(fa.keys())
    chroms = args.chroms if args.chroms else all_chroms
    print(f"Chromosomes: {chroms}")

    chrom_sizes = {c: len(fa[c]) for c in chroms if c in fa}

    # Build sliding-window regions
    regions = []
    for chrom in chroms:
        if chrom not in chrom_sizes:
            print(f"  WARNING: {chrom} not in fasta, skipping")
            continue
        clen = chrom_sizes[chrom]
        for s in range(0, clen - window_size + 1, STRIDE):
            regions.append((chrom, s, s + window_size))
        # last partial window
        if clen > window_size and (clen - window_size) % STRIDE != 0:
            s = clen - window_size
            if not regions or regions[-1] != (chrom, s, s + window_size):
                regions.append((chrom, s, s + window_size))
    print(f"{len(regions):,} windows to predict")

    # Collect per-bin records
    TRACKS = ["G1", "ES", "MS", "LS", "WRT"]
    records = {t: [] for t in TRACKS}

    batch_oh, batch_meta = [], []

    def flush():
        if not batch_oh:
            return
        x = torch.tensor(np.stack(batch_oh), dtype=torch.float32).to(device)
        if not args.no_rc:
            log1p_pred = rc_average(model, x, head="rice").cpu().numpy()
        else:
            with torch.no_grad():
                log1p_pred = model(x, head="rice")["rt_signals"].cpu().numpy()

        pred = np.expm1(np.clip(log1p_pred, 0, None))  # [B, 896, 4] CPM scale
        eps = 1e-6
        g1, es, ms, ls = pred[..., 0], pred[..., 1], pred[..., 2], pred[..., 3]
        wrt = (0.5 * ms + ls) / (es + ms + ls + eps)

        for i, (chrom, win_start, _) in enumerate(batch_meta):
            out_start = win_start + CROP_BINS * BIN_SIZE
            for k in range(OUT_BINS):
                bs = out_start + k * BIN_SIZE
                be = bs + BIN_SIZE
                records["G1"].append((chrom, bs, be, float(g1[i, k])))
                records["ES"].append((chrom, bs, be, float(es[i, k])))
                records["MS"].append((chrom, bs, be, float(ms[i, k])))
                records["LS"].append((chrom, bs, be, float(ls[i, k])))
                records["WRT"].append((chrom, bs, be, float(wrt[i, k])))

        batch_oh.clear()
        batch_meta.clear()

    model.eval()
    for idx, (chrom, ws, we) in enumerate(regions):
        oh = fetch_one_hot(genome, chrom, ws, we, window_size)
        batch_oh.append(oh)
        batch_meta.append((chrom, ws, we))
        if len(batch_oh) == args.batch_size:
            flush()
        if (idx + 1) % 100 == 0:
            print(f"  {idx+1:,}/{len(regions):,} windows ...", flush=True)
    flush()

    # Write BigWigs
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    header = [(c, chrom_sizes[c]) for c in chroms if c in chrom_sizes]

    for track in TRACKS:
        out_path = out_dir / f"{track}.bw"
        entries = sorted(records[track], key=lambda x: (chroms.index(x[0]), x[1]))
        bw = pyBigWig.open(str(out_path), "w")
        bw.addHeader(header)
        for chrom, start, end, val in entries:
            bw.addEntries([chrom], [start], ends=[end], values=[val])
        bw.close()
        print(f"  Written {out_path}")

    print(f"\nDone. {len(TRACKS)} BigWig tracks → {out_dir}")


def main():
    parser = argparse.ArgumentParser(description="Zero-shot cross-species prediction")
    parser.add_argument("--checkpoint", required=True, help="Path to best_model.pt")
    parser.add_argument("--config",     required=True, help="Path to transformer_wt.yaml")
    parser.add_argument("--fasta",      required=True, help="Target species genome fasta")
    parser.add_argument("--species",    required=True, help="Species name (for output labeling)")
    parser.add_argument("--chroms",     nargs="+",     help="Chromosomes to predict (default: all)")
    parser.add_argument("--out_dir",    required=True, help="Output directory for BigWig files")
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--no_rc",      action="store_true", help="Disable RC averaging")
    args = parser.parse_args()
    predict_genome(args)


if __name__ == "__main__":
    main()
