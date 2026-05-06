import random
import numpy as np
import torch
import yaml
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from torch.utils.data import Dataset, Sampler

from .data_utils import GenomeSequence, get_window_coords, load_labels, reverse_complement, RT_CLASS_MAP, one_hot_encode


# ── manifest ──────────────────────────────────────────────────────────────────
@dataclass
class SpeciesConfig:
    name: str
    fasta: str
    count_tsv: str
    gff3: str
    train_chroms: list[str]
    val_chroms: list[str]
    test_chroms: list[str]
    species_id: int


def load_manifest(manifest_yaml: str | Path) -> list[SpeciesConfig]:
    with open(manifest_yaml) as f:
        cfg = yaml.safe_load(f)
    return [
        SpeciesConfig(
            name=sp["name"], fasta=sp["fasta"], count_tsv=sp["count_tsv"],
            gff3=sp["gff3"],
            train_chroms=sp["train_chroms"], val_chroms=sp["val_chroms"],
            test_chroms=sp["test_chroms"], species_id=i,
        )
        for i, sp in enumerate(cfg["species"])
    ]


# ── dataset ───────────────────────────────────────────────────────────────────
class RepliSeqDataset(Dataset):
    def __init__(
        self,
        species_configs: list[SpeciesConfig],
        split: str,                  # "train" | "val" | "test"
        window_size: int = 32768,
        rc_prob: float = 0.5,        # 0.0 disables RC augmentation
        shift_max: int = 0,          # max bp shift each side; 0 disables
    ):
        self.window_size = window_size
        self.rc_prob = rc_prob if split == "train" else 0.0
        self.shift_max = shift_max if split == "train" else 0
        self.samples: list[dict] = []
        self.genomes: dict[str, GenomeSequence] = {}

        for sp in species_configs:
            self.genomes[sp.name] = GenomeSequence(sp.fasta)
            df = load_labels(sp.count_tsv, sp.name, sp.gff3)
            chroms = getattr(sp, f"{split}_chroms")
            df = df[df["chrom"].isin(chroms)].reset_index(drop=True)

            for _, row in df.iterrows():
                cs = self.genomes[sp.name].chrom_size(row["chrom"])
                _, ws, we = get_window_coords(
                    row["chrom"], row["start"], row["end"],
                    window_size=window_size, chrom_size=cs,
                )
                self.samples.append({
                    "species": sp.name, "species_id": sp.species_id,
                    "chrom": row["chrom"], "win_start": ws, "win_end": we,
                    "ES_log1p": float(row["ES_log1p"]),
                    "MS_log1p": float(row["MS_log1p"]),
                    "LS_log1p": float(row["LS_log1p"]),
                    "G1_log1p": float(row["G1_log1p"]),
                    "WRT": float(row["WRT"]),
                    "rt_class": RT_CLASS_MAP[row["RT_class"]],
                })

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> dict:
        s = self.samples[idx]
        shift = random.randint(-self.shift_max, self.shift_max) if self.shift_max > 0 else 0
        seq = self.genomes[s["species"]].fetch(
            s["chrom"], s["win_start"] + shift, s["win_end"] + shift, self.window_size
        )
        if self.rc_prob > 0 and random.random() < self.rc_prob:
            seq = reverse_complement(seq)
        return {
            "one_hot": torch.tensor(one_hot_encode(seq), dtype=torch.float32),  # [4, L]
            "species_id": torch.tensor(s["species_id"], dtype=torch.long),
            "phase_labels": torch.tensor([s["ES_log1p"], s["MS_log1p"], s["LS_log1p"], s["G1_log1p"]], dtype=torch.float32),
            "rt_class": torch.tensor(s["rt_class"], dtype=torch.long),
            "wrt": torch.tensor(s["WRT"], dtype=torch.float32),
        }


# ── sampler ───────────────────────────────────────────────────────────────────
class SpeciesBalancedSampler(Sampler):
    """Equal samples per species per batch; small species are over-sampled."""

    def __init__(self, dataset: RepliSeqDataset, batch_size: int):
        species_idx: dict[str, list[int]] = defaultdict(list)
        for i, s in enumerate(dataset.samples):
            species_idx[s["species"]].append(i)
        self.names = list(species_idx.keys())
        self.idx = {k: np.array(v) for k, v in species_idx.items()}
        self.sps = max(1, batch_size // len(self.names))
        self.n_batches = max(len(v) for v in self.idx.values()) // self.sps

    def __iter__(self):
        shuffled = {k: np.random.permutation(v) for k, v in self.idx.items()}
        ptr = {k: 0 for k in self.names}
        for _ in range(self.n_batches):
            batch = []
            for sp in self.names:
                arr, p, n = shuffled[sp], ptr[sp], self.sps
                if p + n > len(arr):
                    extra = (p + n) - len(arr)
                    chunk = np.concatenate([arr[p:], arr[:extra]])
                    ptr[sp] = extra
                else:
                    chunk = arr[p: p + n]
                    ptr[sp] = p + n
                batch.extend(chunk.tolist())
            yield batch

    def __len__(self) -> int:
        return self.n_batches
