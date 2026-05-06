"""
Genome-wide inference -> BigWig tracks.

BED mode:
  python -m src.infer.predict \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --bed regions.bed --output-dir output/predict/

Chrom mode:
  python -m src.infer.predict \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --chrom chr01 --stride 1000 --output-dir output/predict/
"""
import argparse
import numpy as np
import torch
import pyBigWig
from pathlib import Path
from pyfaidx import Fasta

from ._utils import load_model, parse_bed, fetch_one_hot, rc_average
from ..data.data_utils import GenomeSequence
from ..eval import compute_wrt_from_phase_pred

TRACKS = ["ES", "MS", "LS", "G1", "WRT"]


def _write_bigwigs(out_dir: Path, chrom_sizes: dict, records: dict):
    out_dir.mkdir(parents=True, exist_ok=True)
    for track, entries in records.items():
        entries.sort(key=lambda x: (x[0], x[1]))
        bw = pyBigWig.open(str(out_dir / f"{track}.bw"), "w")
        bw.addHeader(list(chrom_sizes.items()))
        for chrom, start, end, val in entries:
            bw.addEntries([chrom], [start], ends=[end], values=[float(val)])
        bw.close()


def run(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, cfg = load_model(args.checkpoint, args.config, device)
    window_size = cfg["data"]["input_window_length"]
    genome = GenomeSequence(args.fasta)
    fa = Fasta(args.fasta)
    chrom_sizes = {c: len(fa[c]) for c in fa.keys()}

    records = {t: [] for t in TRACKS}

    if args.bed:
        regions = list(parse_bed(args.bed))
    else:
        chrom = args.chrom
        clen = chrom_sizes[chrom]
        stride = args.stride
        regions = [
            (chrom, s, s + stride, f"{chrom}:{s}-{s+stride}", "+")
            for s in range(0, clen - stride + 1, stride)
        ]

    batch_oh, batch_meta = [], []

    def flush(batch_oh, batch_meta):
        if not batch_oh:
            return
        x = torch.tensor(np.stack(batch_oh), dtype=torch.float32).to(device)
        if args.rc:
            pred = rc_average(model, x).cpu().numpy()
        else:
            with torch.no_grad():
                pred = model(x)["phase_pred"].cpu().numpy()
        wrt = compute_wrt_from_phase_pred(pred)
        for i, (chrom, start, end) in enumerate(batch_meta):
            for j, t in enumerate(["ES", "MS", "LS", "G1"]):
                records[t].append((chrom, start, end, pred[i, j]))
            records["WRT"].append((chrom, start, end, wrt[i]))

    for chrom, start, end, name, strand in regions:
        oh = fetch_one_hot(genome, chrom, start, end, window_size, strand)
        batch_oh.append(oh)
        batch_meta.append((chrom, start, end))
        if len(batch_oh) == args.batch_size:
            flush(batch_oh, batch_meta)
            batch_oh, batch_meta = [], []
    flush(batch_oh, batch_meta)

    _write_bigwigs(Path(args.output_dir), chrom_sizes, records)
    print(f"Written {len(TRACKS)} BigWig tracks to {args.output_dir}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--bed")
    p.add_argument("--chrom")
    p.add_argument("--stride", type=int, default=1000)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--rc", action="store_true")
    p.add_argument("--output-dir", default="output/predict")
    args = p.parse_args()
    assert args.bed or args.chrom, "Provide --bed or --chrom"
    run(args)


if __name__ == "__main__":
    main()
