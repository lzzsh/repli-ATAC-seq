"""
In silico saturation mutagenesis -> HDF5.

  python -m src.infer.ism \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --bed regions.bed --output ism_scores.h5

  python -m src.infer.ism \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --region chr01:1000000-1001000 --output ism_scores.h5
"""
import argparse
import numpy as np
import torch
import h5py
from pathlib import Path

from ._utils import load_model, parse_bed, parse_region, fetch_one_hot
from ..data.data_utils import GenomeSequence


def ism_region(model, one_hot_np: np.ndarray, device: torch.device,
               mut_batch_size: int):
    """
    Returns:
        ref_pred:   [4]             log1p(TPM) predictions for reference
        ism_scores: [L, 4, 4]      mean-normalized (pred_mut - pred_ref), axes: [pos, base, phase]
        scd:        [L, 4]         sqrt(sum(ism_scores^2, axis=1)) per position
    """
    L = one_hot_np.shape[1]
    ref_t = torch.tensor(one_hot_np[None], dtype=torch.float32).to(device)
    with torch.no_grad():
        ref_pred = model(ref_t)["phase_pred"].squeeze(0).cpu().numpy()  # [4]

    ism_scores = np.zeros((L, 4, 4), dtype=np.float32)

    mut_list = []
    for i in range(L):
        ref_base = one_hot_np[:, i].argmax() if one_hot_np[:, i].sum() > 0 else -1
        for b in range(4):
            if b == ref_base:
                continue
            mut = one_hot_np.copy()
            mut[:, i] = 0.0
            mut[b, i] = 1.0
            mut_list.append((i, b, mut))

    for batch_start in range(0, len(mut_list), mut_batch_size):
        batch = mut_list[batch_start: batch_start + mut_batch_size]
        x = torch.tensor(
            np.stack([m[2] for m in batch]), dtype=torch.float32
        ).to(device)
        with torch.no_grad():
            preds = model(x)["phase_pred"].cpu().numpy()  # [B, 4]
        for k, (i, b, _) in enumerate(batch):
            ism_scores[i, b, :] = preds[k] - ref_pred

    ism_scores -= ism_scores.mean(axis=1, keepdims=True)
    scd = np.sqrt((ism_scores ** 2).sum(axis=1))  # [L, 4]
    return ref_pred, ism_scores, scd


def run(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, cfg = load_model(args.checkpoint, args.config, device)
    window_size = cfg["data"]["input_window_length"]
    genome = GenomeSequence(args.fasta)

    if args.bed:
        regions = list(parse_bed(args.bed))
    else:
        chrom, start, end = parse_region(args.region)
        regions = [(chrom, start, end, f"{chrom}:{start}-{end}", "+")]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    n = len(regions)
    L = window_size

    with h5py.File(args.output, "w") as hf:
        hf.create_dataset("seqnames", data=np.array([r[3] for r in regions], dtype="S"))
        hf.create_dataset("chrom",    data=np.array([r[0] for r in regions], dtype="S"))
        hf.create_dataset("start",    data=np.array([r[1] for r in regions], dtype=np.int32))
        hf.create_dataset("end",      data=np.array([r[2] for r in regions], dtype=np.int32))
        hf.create_dataset("strand",   data=np.array([r[4] for r in regions], dtype="S"))
        ds_seqs = hf.create_dataset("seqs",       shape=(n, 4, L),    dtype=np.float32)
        ds_ref  = hf.create_dataset("ref_pred",   shape=(n, 4),       dtype=np.float32)
        ds_ism  = hf.create_dataset("ism_scores", shape=(n, L, 4, 4), dtype=np.float32)
        ds_scd  = hf.create_dataset("scd",        shape=(n, L, 4),    dtype=np.float32)

        for idx, (chrom, start, end, name, strand) in enumerate(regions):
            print(f"[{idx+1}/{n}] {name}")
            oh = fetch_one_hot(genome, chrom, start, end, window_size, strand)
            ref_pred, ism_scores, scd = ism_region(model, oh, device, args.mut_batch_size)
            ds_seqs[idx] = oh
            ds_ref[idx]  = ref_pred
            ds_ism[idx]  = ism_scores
            ds_scd[idx]  = scd

    print(f"Written ISM scores to {args.output}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--bed")
    p.add_argument("--region")
    p.add_argument("--mut-batch-size", type=int, default=64)
    p.add_argument("--output", default="ism_scores.h5")
    args = p.parse_args()
    assert args.bed or args.region, "Provide --bed or --region"
    run(args)


if __name__ == "__main__":
    main()
