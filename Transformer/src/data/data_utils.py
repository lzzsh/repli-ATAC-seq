import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

# ── constants ────────────────────────────────────────────────────────────────
RT_CLASS_MAP = {"ES": 0, "MS": 1, "LS": 2, "NR": 3}   # NR = not replicated
VALID_GFF3_CLASSES = {"ES", "MS", "LS"}
BOUNDARY_CLASSES = {"ESMS", "MSLS", "ESLS", "ESMSLS"}
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
) -> tuple[dict[tuple[str, int, int], str], set[tuple[str, int, int]]]:
    """Parse GFF3 in one pass.

    Returns:
        classes: (chrom, start, end) → Name for ES/MS/LS bins only
        boundary_keys: (chrom, start, end) for boundary bins to be filtered
    """
    classes: dict[tuple[str, int, int], str] = {}
    boundary_keys: set[tuple[str, int, int]] = set()
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, start, end = parts[0], int(float(parts[3])), int(float(parts[4]))
            name = ""
            for attr in parts[8].split(";"):
                if attr.startswith("Name="):
                    name = attr[5:]
                    break
            if name in BOUNDARY_CLASSES:
                boundary_keys.add((chrom, start, end))
            elif name in VALID_GFF3_CLASSES:
                classes[(chrom, start, end)] = name
    return classes, boundary_keys

def load_labels(count_tsv: str | Path, species: str, gff3_path: str | Path) -> pd.DataFrame:
    """
    Input TSV: chrom, start, end, ES_count, MS_count, LS_count, G1_count.
    GFF3: replication phase annotations (ES/MS/LS/boundary/absent).
    Returns DataFrame with TPM, log1p, WRT columns and RT_class in {ES,MS,LS,NR}.
    Boundary bins (ESMS/MSLS/ESLS/ESMSLS) are filtered out.
    """
    df = pd.read_csv(count_tsv, sep="\t")
    for col in ["chrom", "start", "end", "ES_count", "MS_count", "LS_count", "G1_count"]:
        assert col in df.columns, f"Missing column: {col}"

    es = tpm_normalize(df["ES_count"].values)
    ms = tpm_normalize(df["MS_count"].values)
    ls = tpm_normalize(df["LS_count"].values)
    g1 = tpm_normalize(df["G1_count"].values)

    df["ES_tpm"], df["MS_tpm"], df["LS_tpm"], df["G1_tpm"] = es, ms, ls, g1
    df["ES_log1p"] = np.log1p(es).astype(np.float32)
    df["MS_log1p"] = np.log1p(ms).astype(np.float32)
    df["LS_log1p"] = np.log1p(ls).astype(np.float32)
    df["G1_log1p"] = np.log1p(g1).astype(np.float32)
    df["WRT"] = compute_wrt(es, ms, ls)
    df["species"] = species

    gff3_classes, boundary_keys = load_gff3_classes(gff3_path)
    df["RT_class"] = df.apply(
        lambda r: gff3_classes.get((r["chrom"], r["start"], r["end"]), "NR"), axis=1
    )
    is_boundary = df.apply(
        lambda r: (r["chrom"], r["start"], r["end"]) in boundary_keys, axis=1
    )
    return df[~is_boundary].reset_index(drop=True)
