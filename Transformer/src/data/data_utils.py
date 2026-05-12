import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

# ── constants ────────────────────────────────────────────────────────────────
RT_CLASS_MAP = {"ES": 0, "MS": 1, "LS": 2, "NR": 3}
VALID_GFF3_CLASSES = {"ES", "MS", "LS", "NR"}
# IGV display names → internal class names
_GFF3_NAME_MAP = {"Non-replication": "NR", "unknown": "NR"}
_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


# ── sequence utilities ────────────────────────────────────────────────────────
def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def one_hot_encode(seq: str) -> np.ndarray:
    """Encode DNA string to one-hot array of shape [4, L].

    Channels: A=0, C=1, G=2, T=3. Non-ACGT bases (N, IUPAC) → all zeros.
    """
    L = len(seq)
    out = np.zeros((4, L), dtype=np.float32)
    for i, base in enumerate(seq.upper()):
        idx = _BASE_IDX.get(base)
        if idx is not None:
            out[idx, i] = 1.0
    return out


class GenomeSequence:
    def __init__(self, fasta_path: str | Path):
        self.fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=False)

    def fetch(self, chrom: str, start: int, end: int, window_size: int = 32768) -> str:
        """Fetch [start, end), pad with N at boundaries, soft-mask → uppercase."""
        chrom_len = len(self.fasta[chrom])
        pad_left = max(0, -start)
        pad_right = max(0, end - chrom_len)
        seq = str(self.fasta[chrom][max(0, start):min(chrom_len, end)]).upper()
        seq = "N" * pad_left + seq + "N" * pad_right
        seq = seq + "N" * max(0, window_size - len(seq))
        return seq[:window_size]

    def chrom_size(self, chrom: str) -> int:
        return len(self.fasta[chrom])


# ── bin / window utilities ────────────────────────────────────────────────────
def get_window_coords(
    chrom: str, bin_start: int, bin_end: int,
    window_size: int = 32768, chrom_size: int | None = None,
) -> tuple[str, int, int]:
    """Center an 8 kb window on the bin; clamp to chromosome bounds."""
    center = (bin_start + bin_end) // 2
    half = window_size // 2
    ws, we = center - half, center + half
    if chrom_size is not None:
        ws, we = max(0, ws), min(chrom_size, we)
    return chrom, ws, we


# ── normalization & labels ────────────────────────────────────────────────────
def tpm_normalize(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts, dtype=np.float32)
    return (counts / total * 1e6).astype(np.float32)


def compute_wrt(es: np.ndarray, ms: np.ndarray, ls: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return (0.5 * ms + ls) / (es + ms + ls + eps)


def load_gff3_classes(
    gff3_path: str | Path,
) -> dict[tuple[str, int, int], str]:
    """Parse GFF3, return (chrom, start, end) → RT class name for all valid classes."""
    classes: dict[tuple[str, int, int], str] = {}
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, start, end = parts[0], int(float(parts[3])) - 1, int(float(parts[4]))
            name = ""
            for attr in parts[8].split(";"):
                if attr.startswith("Name="):
                    name = attr[5:].strip()
                    break
            name = _GFF3_NAME_MAP.get(name, name)
            if name in VALID_GFF3_CLASSES:
                classes[(chrom, start, end)] = name
    return classes

def load_labels(gff3_path: str | Path, species: str) -> pd.DataFrame:
    """
    Parse GFF3 replication phase annotations into a DataFrame with columns:
    chrom, start, end, RT_class, species.
    """
    gff3_classes = load_gff3_classes(gff3_path)
    rows = [
        {"chrom": chrom, "start": start, "end": end, "RT_class": cls, "species": species}
        for (chrom, start, end), cls in gff3_classes.items()
    ]
    return pd.DataFrame(rows).reset_index(drop=True)


IGNORE_LABEL = -1   # bins with no GFF3 annotation; excluded from loss and metrics


def load_labels_indexed(df: pd.DataFrame):
    """
    Returns query(chrom, genomic_start, n_bins, model_bin_size) -> np.ndarray [n_bins] int64
    Each model bin maps to the label bin containing its center coordinate.
    Bins with no annotation are set to IGNORE_LABEL (-1), not NR.
    """
    grouped = {}
    for chrom, grp in df.groupby("chrom"):
        label_bin_size = int(grp["end"].iloc[0] - grp["start"].iloc[0])
        keys = grp["start"].values.astype(int)
        vals = grp["RT_class"].map(RT_CLASS_MAP).fillna(IGNORE_LABEL).values.astype(np.int64)
        idx_map = dict(zip(keys, vals))
        grouped[chrom] = (label_bin_size, idx_map)

    def query(chrom: str, genomic_start: int, n_bins: int, model_bin_size: int) -> np.ndarray:
        out = np.full(n_bins, IGNORE_LABEL, dtype=np.int64)
        if chrom not in grouped:
            return out
        label_bin_size, idx_map = grouped[chrom]
        for i in range(n_bins):
            center = genomic_start + i * model_bin_size + model_bin_size // 2
            label_key = (center // label_bin_size) * label_bin_size
            val = idx_map.get(label_key)
            if val is not None:
                out[i] = val
        return out

    return query
