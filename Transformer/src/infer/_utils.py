import numpy as np
import torch
import yaml
from pathlib import Path

from ..data.data_utils import GenomeSequence, one_hot_encode, reverse_complement
from ..data.dataset import load_manifest
from ..models.model import RepliformerModel
from ..models.config_model import RepliformerConfig


def load_model(checkpoint: str, config: str, device: torch.device):
    with open(config) as f:
        cfg = yaml.safe_load(f)
    manifest_path = cfg["data"]["manifest"]
    if not Path(manifest_path).is_absolute():
        # try progressively higher parent directories until manifest is found
        config_path = Path(config).resolve()
        for parent in config_path.parents:
            candidate = parent / manifest_path
            if candidate.exists():
                manifest_path = str(candidate)
                break
    species_configs = load_manifest(manifest_path)
    ckpt = torch.load(checkpoint, map_location="cpu")
    # use model config saved in checkpoint if available, else fall back to yaml
    ckpt_cfg = ckpt.get("cfg", {})
    model_cfg_dict = ckpt_cfg.get("model", cfg.get("model", {}))
    model = RepliformerModel(
        species_configs=species_configs,
        model_cfg=RepliformerConfig(**model_cfg_dict),
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval(), ckpt_cfg or cfg


def parse_bed(bed_path: str):
    """Yield (chrom, start, end, name, strand) from BED file."""
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            strand = parts[5] if len(parts) > 5 else "+"
            yield chrom, start, end, name, strand


def parse_region(region_str: str):
    """Parse 'chr:start-end' string -> (chrom, start, end)."""
    chrom, rest = region_str.split(":")
    start, end = rest.split("-")
    return chrom, int(start), int(end)


def fetch_one_hot(genome: GenomeSequence, chrom: str, start: int, end: int,
                  window_size: int, strand: str = "+") -> np.ndarray:
    """Center window on [start,end), fetch sequence, one-hot encode. Shape [4, L]."""
    center = (start + end) // 2
    half = window_size // 2
    ws = center - half
    we = center + half
    chrom_size = genome.chrom_size(chrom)
    ws, we = max(0, ws), min(chrom_size, we)
    seq = genome.fetch(chrom, ws, we, window_size)
    if strand == "-":
        seq = reverse_complement(seq)
    return one_hot_encode(seq)  # [4, L]


def rc_average(model: RepliformerModel, batch: torch.Tensor, head: str) -> torch.Tensor:
    """Average forward and reverse-complement predictions. batch: [B, 4, L]."""
    with torch.no_grad():
        fwd = model(batch, head=head)["rt_signals"]                        # [B, T, 4]
        rc  = torch.flip(batch, dims=[-1])[:, [3, 2, 1, 0], :]
        rev = torch.flip(model(rc, head=head)["rt_signals"], dims=[1])     # flip bin dim
    return (fwd + rev) / 2
