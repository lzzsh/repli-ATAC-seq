"""
Gradient saliency -> HDF5.

  python -m src.infer.saliency \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --bed regions.bed --output saliency_scores.h5

  python -m src.infer.saliency \
    --checkpoint best_model.pt --config src/configs/transformer_wt.yaml \
    --fasta genome.fa --region chr01:1000000-1001000 \
    --smooth-samples 50 --smooth-noise 0.1 --output saliency_scores.h5
"""
import argparse
import numpy as np
import torch
import h5py
from pathlib import Path

from ._utils import load_model, parse_bed, parse_region, fetch_one_hot
from ..data.data_utils import GenomeSequence

N_PHASES = 4


def saliency_region(model, one_hot_np: np.ndarray, device: torch.device,
                    smooth_samples: int, smooth_noise: float):
    """
    Returns:
        ref_pred:    [4]               log1p(TPM) predictions
        saliency:    [N_PHASES, L, 4]  grad x input  (axes: phase, position, base)
        hypothetical:[N_PHASES, L, 4]  grad only     (axes: phase, position, base)
    """
    L = one_hot_np.shape[1]
    base_t = torch.tensor(one_hot_np[None], dtype=torch.float32).to(device)

    with torch.no_grad():
        ref_pred = model(base_t)["phase_pred"].squeeze(0).cpu().numpy()

    def _compute_grad(x_t: torch.Tensor) -> np.ndarray:
        grads = []
        for ph in range(N_PHASES):
            model.zero_grad()
            if x_t.grad is not None:
                x_t.grad.zero_()
            pred = model(x_t)["phase_pred"]
            pred[0, ph].backward(retain_graph=(ph < N_PHASES - 1))
            grads.append(x_t.grad.detach().squeeze(0).cpu().numpy())  # [4, L]
        return np.stack(grads)  # [N_PHASES, 4, L]

    if smooth_samples > 0:
        acc = np.zeros((N_PHASES, 4, L), dtype=np.float32)
        for _ in range(smooth_samples):
            noise = torch.randn_like(base_t) * smooth_noise
            x_t = (base_t + noise).requires_grad_(True)
            acc += _compute_grad(x_t)
        grad = acc / smooth_samples
    else:
        x_t = base_t.requires_grad_(True)
        grad = _compute_grad(x_t)  # [N_PHASES, 4, L]

    # transpose [N_PHASES, 4, L] -> [N_PHASES, L, 4]
    grad = grad.transpose(0, 2, 1)
    oh_t = one_hot_np.T[None]        # [1, L, 4]
    saliency     = grad * oh_t       # [N_PHASES, L, 4]
    hypothetical = grad.copy()       # [N_PHASES, L, 4]

    return ref_pred, saliency, hypothetical


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
        ds_seqs = hf.create_dataset("seqs",         shape=(n, 4, L),             dtype=np.float32)
        ds_ref  = hf.create_dataset("ref_pred",     shape=(n, 4),                dtype=np.float32)
        ds_sal  = hf.create_dataset("saliency",     shape=(n, N_PHASES, L, 4),   dtype=np.float32)
        ds_hyp  = hf.create_dataset("hypothetical", shape=(n, N_PHASES, L, 4),   dtype=np.float32)

        for idx, (chrom, start, end, name, strand) in enumerate(regions):
            print(f"[{idx+1}/{n}] {name}")
            oh = fetch_one_hot(genome, chrom, start, end, window_size, strand)
            ref_pred, saliency, hypothetical = saliency_region(
                model, oh, device, args.smooth_samples, args.smooth_noise
            )
            ds_seqs[idx] = oh
            ds_ref[idx]  = ref_pred
            ds_sal[idx]  = saliency
            ds_hyp[idx]  = hypothetical

    print(f"Written saliency scores to {args.output}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--bed")
    p.add_argument("--region")
    p.add_argument("--smooth-samples", type=int, default=0)
    p.add_argument("--smooth-noise",   type=float, default=0.1)
    p.add_argument("--output", default="saliency_scores.h5")
    args = p.parse_args()
    assert args.bed or args.region, "Provide --bed or --region"
    run(args)


if __name__ == "__main__":
    main()
