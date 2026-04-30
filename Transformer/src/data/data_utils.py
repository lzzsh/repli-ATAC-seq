import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

# ── constants ────────────────────────────────────────────────────────────────
RT_CLASS_MAP = {"E": 0, "M": 1, "L": 2}
VALID_RT_CLASSES = {"E", "M", "L"}
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

    def fetch(self, chrom: str, start: int, end: int, window_size: int = 8192) -> str:
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
    window_size: int = 8192, chrom_size: int | None = None,
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


def assign_rt_class(wrt: np.ndarray, low: float = 1 / 3, high: float = 2 / 3) -> np.ndarray:
    cls = np.full(len(wrt), "M", dtype=object)
    cls[wrt < low] = "E"
    cls[wrt > high] = "L"
    return cls


def load_labels(count_tsv: str | Path, species: str) -> pd.DataFrame:
    """
    Input TSV: chrom, start, end, ES_count, MS_count, LS_count.
    Returns DataFrame with TPM, log1p, WRT, RT_class columns.
    Only rows with RT_class in {E, M, L} are kept (EM/ML excluded).
    """
    df = pd.read_csv(count_tsv, sep="\t")
    for col in ["chrom", "start", "end", "ES_count", "MS_count", "LS_count"]:
        assert col in df.columns, f"Missing column: {col}"

    es = tpm_normalize(df["ES_count"].values)
    ms = tpm_normalize(df["MS_count"].values)
    ls = tpm_normalize(df["LS_count"].values)

    df["ES_tpm"], df["MS_tpm"], df["LS_tpm"] = es, ms, ls
    df["ES_log1p"] = np.log1p(es).astype(np.float32)
    df["MS_log1p"] = np.log1p(ms).astype(np.float32)
    df["LS_log1p"] = np.log1p(ls).astype(np.float32)
    df["WRT"] = compute_wrt(es, ms, ls)
    df["RT_class"] = assign_rt_class(df["WRT"].values)
    df["species"] = species
    return df[df["RT_class"].isin(VALID_RT_CLASSES)].reset_index(drop=True)
