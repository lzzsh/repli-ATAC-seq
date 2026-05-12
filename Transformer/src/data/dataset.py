import random
import numpy as np
import torch
import yaml
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from torch.utils.data import Dataset, Sampler

from .data_utils import GenomeSequence, get_window_coords, load_labels, load_labels_indexed, reverse_complement, RT_CLASS_MAP, one_hot_encode


# ── manifest ──────────────────────────────────────────────────────────────────
@dataclass
class SpeciesConfig:
    name: str
    fasta: str
    gff3: str
    train_chroms: list[str]
    val_chroms: list[str]
    test_chroms: list[str]
    species_id: int


def load_manifest(manifest_yaml: str | Path) -> list[SpeciesConfig]:
    with open(manifest_yaml) as f:
        cfg = yaml.safe_load(f)
    ids = [sp["species_id"] for sp in cfg["species"]]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate species_id values in manifest: {ids}")
    return [
        SpeciesConfig(
            name=sp["name"], fasta=sp["fasta"],
            gff3=sp["gff3"],
            train_chroms=sp["train_chroms"], val_chroms=sp["val_chroms"],
            test_chroms=sp["test_chroms"], species_id=sp["species_id"],
        )
        for sp in cfg["species"]
    ]


# ── dataset ───────────────────────────────────────────────────────────────────
class RepliSeqDataset(Dataset):
    _BIN_SIZE    = 128      # bp per trunk output bin
    _OUT_BIN     = 128      # bp per model output bin (no pooling)
    _OUT_BINS    = 896      # 196608/128 - 640 = 896 (crop 320 each side)
    _CROP_BINS   = 320      # trunk crops 320 bins each side at 128bp (= 40960 bp)
    _STRIDE      = 114688   # bp between window starts (= output region, zero overlap, matches Enformer)

    def __init__(
        self,
        species_configs: list[SpeciesConfig],
        split: str,
        window_size: int = 32768,
        rc_prob: float = 0.5,
        shift_max: int = 0,
    ):
        self.window_size = window_size
        self.rc_prob = rc_prob if split == "train" else 0.0
        self.shift_max = shift_max if split == "train" else 0
        self.samples: list[dict] = []
        self.genomes: dict[str, GenomeSequence] = {}
        self._label_queries: dict[str, object] = {}

        for sp in species_configs:
            self.genomes[sp.name] = GenomeSequence(sp.fasta)
            df = load_labels(sp.gff3, sp.name)
            chroms = getattr(sp, f"{split}_chroms")
            df = df[df["chrom"].isin(chroms)].reset_index(drop=True)
            self._label_queries[sp.name] = load_labels_indexed(df)

            annotated = set(zip(df["chrom"], df["start"]))
            for chrom in chroms:
                chrom_size = self.genomes[sp.name].chrom_size(chrom)
                for win_start in range(0, chrom_size - window_size + 1, self._STRIDE):
                    # first output bin starts after the 320-bin trunk crop
                    out_offset = win_start + self._CROP_BINS * self._BIN_SIZE
                    has_label = any(
                        (chrom, ((out_offset + i * self._OUT_BIN + self._OUT_BIN // 2)
                                 // self._OUT_BIN) * self._OUT_BIN)
                        in annotated
                        for i in range(self._OUT_BINS)
                    )
                    if has_label:
                        self.samples.append({
                            "species": sp.name,
                            "species_id": sp.species_id,
                            "chrom": chrom,
                            "win_start": win_start,
                        })

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> dict:
        s = self.samples[idx]
        shift = random.randint(-self.shift_max, self.shift_max) if self.shift_max > 0 else 0
        win_start = s["win_start"] + shift
        seq = self.genomes[s["species"]].fetch(
            s["chrom"], win_start, win_start + self.window_size, self.window_size
        )
        rc = self.rc_prob > 0 and random.random() < self.rc_prob
        if rc:
            seq = reverse_complement(seq)

        out_offset = win_start + self._CROP_BINS * self._BIN_SIZE
        labels = self._label_queries[s["species"]](
            s["chrom"], out_offset, self._OUT_BINS, self._OUT_BIN
        )  # [896] int64
        if rc:
            labels = labels[::-1].copy()

        return {
            "one_hot": torch.tensor(one_hot_encode(seq), dtype=torch.float32),
            "species_id": torch.tensor(s["species_id"], dtype=torch.long),
            "rt_labels": torch.tensor(labels, dtype=torch.long),  # [896]
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
