import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta
from typing import Optional

# ── constants ────────────────────────────────────────────────────────────────
SIGNAL_CHANNELS = ["G1", "ES", "MS", "LS"]
N_SIGNALS = len(SIGNAL_CHANNELS)

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


# ── sequence utilities ────────────────────────────────────────────────────────
def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def one_hot_encode(seq: str) -> np.ndarray:
    """Encode DNA string to one-hot array of shape [4, L].

    Channels: A=0, C=1, G=2, T=3. Non-ACGT bases → all zeros.
    """
    seq_bytes = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    out = np.zeros((4, len(seq_bytes)), dtype=np.float32)
    for i, b in enumerate(b"ACGT"):
        out[i] = (seq_bytes == b)
    return out


class GenomeSequence:
    def __init__(self, fasta_path: str | Path, open_regions_gff3: Optional[str | Path] = None):
        self.fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=False)
        self._open_mask: dict[str, np.ndarray] = {}
        if open_regions_gff3 is not None:
            self._build_open_mask(open_regions_gff3)

    def _build_open_mask(self, gff3_path: str | Path) -> None:
        # GFF3 cols: seqname(0), source(1), feature(2), start(3), end(4), ...
        # GFF3 is 1-based inclusive → convert to 0-based half-open [start-1, end)
        df = pd.read_csv(
            gff3_path, sep="\t", header=None, comment="#",
            usecols=[0, 3, 4], names=["chrom", "start", "end"],
        )
        for chrom, grp in df.groupby("chrom"):
            if chrom not in self.fasta:
                continue
            chrom_len = len(self.fasta[chrom])
            mask = np.zeros(chrom_len, dtype=bool)
            starts = (grp["start"].values - 1).clip(0, chrom_len)   # 1-based → 0-based
            ends   = grp["end"].values.clip(0, chrom_len)
            for s, e in zip(starts, ends):
                mask[s:e] = True
            self._open_mask[chrom] = mask

    def fetch(self, chrom: str, start: int, end: int, window_size: int = 32768) -> str:
        """Fetch [start, end), pad with N at boundaries, soft-mask → uppercase."""
        chrom_len = len(self.fasta[chrom])
        pad_left = max(0, -start)
        pad_right = max(0, end - chrom_len)
        seq = str(self.fasta[chrom][max(0, start):min(chrom_len, end)]).upper()
        seq = "N" * pad_left + seq + "N" * pad_right
        seq = seq + "N" * max(0, window_size - len(seq))
        seq = seq[:window_size]

        if chrom in self._open_mask:
            seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            actual_start = max(0, start)
            actual_end   = min(chrom_len, end)
            genome_pos   = np.arange(actual_start, actual_end)
            win_pos      = np.arange(pad_left, pad_left + len(genome_pos))
            closed       = ~self._open_mask[chrom][genome_pos]
            seq_arr[win_pos[closed]] = ord("N")
            seq = seq_arr.tobytes().decode("ascii")

        return seq

    def chrom_size(self, chrom: str) -> int:
        return len(self.fasta[chrom])


# ── signal label loading ──────────────────────────────────────────────────────
def load_signals(tsv_path: str | Path, species: str) -> pd.DataFrame:
    """
    Parse TSV replication signal file.
    Required columns: chrom, start, end, G1, ES, MS, LS (CPM values).
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype={"chrom": str})
    required = {"chrom", "start", "end"} | set(SIGNAL_CHANNELS)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing columns: {missing}")
    df = df[["chrom", "start", "end"] + SIGNAL_CHANNELS].copy()
    df["species"] = species
    return df.reset_index(drop=True)


def load_signals_indexed(df: pd.DataFrame):
    """
    Returns query(chrom, genomic_start, n_bins, model_bin_size) -> float32 [n_bins, 4]
    Each model bin maps to the signal bin containing its center coordinate.
    Unannotated bins are NaN (excluded from loss via mask).
    """
    grouped = {}
    for chrom, grp in df.groupby("chrom"):
        starts = grp["start"].values.astype(int)
        ends   = grp["end"].values.astype(int)
        vals   = grp[SIGNAL_CHANNELS].values.astype(np.float32)
        order  = np.argsort(starts)
        grouped[chrom] = (starts[order], ends[order], vals[order])

    def query(chrom: str, genomic_start: int, n_bins: int, model_bin_size: int) -> np.ndarray:
        out = np.full((n_bins, N_SIGNALS), np.nan, dtype=np.float32)
        if chrom not in grouped:
            return out
        starts, ends, vals = grouped[chrom]
        centers = genomic_start + np.arange(n_bins) * model_bin_size + model_bin_size // 2
        idxs = np.searchsorted(starts, centers, side="right") - 1
        valid = (idxs >= 0) & (idxs < len(starts))
        valid[valid] &= (starts[idxs[valid]] <= centers[valid]) & (centers[valid] < ends[idxs[valid]])
        out[valid] = vals[idxs[valid]]
        return out

    return query
