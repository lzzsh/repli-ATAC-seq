import random
import numpy as np
import torch
import yaml
from dataclasses import dataclass
from pathlib import Path
from torch.utils.data import Dataset

from .data_utils import GenomeSequence, load_signals, load_signals_indexed, reverse_complement, one_hot_encode


# ── manifest ──────────────────────────────────────────────────────────────────
@dataclass
class SpeciesConfig:
    name: str
    fasta: str
    tsv: str
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
            tsv=sp["tsv"],
            train_chroms=sp["train_chroms"], val_chroms=sp["val_chroms"],
            test_chroms=sp["test_chroms"], species_id=sp["species_id"],
        )
        for sp in cfg["species"]
    ]


# ── dataset ───────────────────────────────────────────────────────────────────
class RepliSeqDataset(Dataset):
    _BIN_SIZE    = 128
    _OUT_BIN     = 128
    _OUT_BINS    = 896
    _CROP_BINS   = 320

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
        self._stride = window_size - 2 * self._CROP_BINS * self._BIN_SIZE
        self.samples: list[dict] = []
        self.genomes: dict[str, GenomeSequence] = {}
        self._signal_queries: dict[str, object] = {}

        for sp in species_configs:
            self.genomes[sp.name] = GenomeSequence(sp.fasta)
            df = load_signals(sp.tsv, sp.name)
            chroms = getattr(sp, f"{split}_chroms")
            df = df[df["chrom"].isin(chroms)].reset_index(drop=True)
            self._signal_queries[sp.name] = load_signals_indexed(df)

            annotated_chroms = set(df["chrom"].unique())
            for chrom in chroms:
                if chrom not in annotated_chroms:
                    continue
                chrom_size = self.genomes[sp.name].chrom_size(chrom)
                for win_start in range(0, chrom_size - window_size + 1, self._stride):
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
        chrom_size = self.genomes[s["species"]].chrom_size(s["chrom"])
        win_start = max(0, min(s["win_start"] + shift, chrom_size - self.window_size))
        seq = self.genomes[s["species"]].fetch(
            s["chrom"], win_start, win_start + self.window_size, self.window_size
        )
        rc = self.rc_prob > 0 and random.random() < self.rc_prob
        if rc:
            seq = reverse_complement(seq)

        out_offset = win_start + self._CROP_BINS * self._BIN_SIZE
        signals = self._signal_queries[s["species"]](
            s["chrom"], out_offset, self._OUT_BINS, self._OUT_BIN
        )  # [896, 4] float32, NaN where unannotated

        # apply log1p to valid (non-NaN) bins; NaN bins stay NaN
        nan_mask = np.isnan(signals)
        signals = np.log1p(np.where(nan_mask, 0.0, signals))
        signals[nan_mask] = np.nan

        if rc:
            signals = signals[::-1].copy()

        return {
            "one_hot": torch.tensor(one_hot_encode(seq), dtype=torch.float32),
            "species_id": torch.tensor(s["species_id"], dtype=torch.long),
            "rt_signals": torch.tensor(signals, dtype=torch.float32),  # [896, 4]
        }
