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
    # trunk outputs 224 bins at 128bp resolution after cropping 16 each side.
    # Center 128 bins (48:176) = 16384bp centered on the 32kb window are used for loss.
    # Stride = 16384bp so adjacent samples' center windows do not overlap.
    _BIN_SIZE     = 128    # bp per output bin
    _CENTER_START = 48     # first of the 128 center bins in cropped space
    _CENTER_BINS  = 128    # bins used for loss / labels  (128 × 128bp = 16384bp)
    _STRIDE       = _CENTER_BINS * _BIN_SIZE  # 16384bp

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
            df = load_labels(sp.count_tsv, sp.name, sp.gff3)
            chroms = getattr(sp, f"{split}_chroms")
            df = df[df["chrom"].isin(chroms)].reset_index(drop=True)
            self._label_queries[sp.name] = load_labels_indexed(df)

            # stride = 16384bp (= _CENTER_BINS × _BIN_SIZE), non-overlapping center windows.
            # Only keep windows whose center 128-bin region has at least one annotated bin.
            annotated = set(zip(df["chrom"], df["start"]))
            for chrom in chroms:
                chrom_size = self.genomes[sp.name].chrom_size(chrom)
                for win_start in range(0, chrom_size - window_size + 1, self._STRIDE):
                    center_offset = win_start + (16 + self._CENTER_START) * self._BIN_SIZE
                    # check if any of the 128 center bins has a label
                    has_label = any(
                        (chrom, (center_offset + i * self._BIN_SIZE) // self._BIN_SIZE * self._BIN_SIZE)
                        in annotated
                        for i in range(self._CENTER_BINS)
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

        # center 128 output bins start at: win_start + (16 crop bins + 48) × 128bp
        center_offset = win_start + (16 + self._CENTER_START) * self._BIN_SIZE
        labels = self._label_queries[s["species"]](
            s["chrom"], center_offset, self._CENTER_BINS, self._BIN_SIZE
        )  # [128, 4]
        if rc:
            labels = labels[::-1].copy()

        return {
            "one_hot": torch.tensor(one_hot_encode(seq), dtype=torch.float32),
            "species_id": torch.tensor(s["species_id"], dtype=torch.long),
            "phase_labels": torch.tensor(labels, dtype=torch.float32),  # [128, 4]
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
