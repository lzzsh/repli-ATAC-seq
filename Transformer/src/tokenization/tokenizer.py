from itertools import product

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _rc_kmer(kmer: str) -> str:
    return "".join(_COMP[b] for b in reversed(kmer))


def build_canonical_vocab(k: int = 6) -> dict[str, int]:
    """
    Special tokens: [PAD]=0, [UNK]=1, [MASK]=2, [CLS]=3.
    Canonical k-mer = lexicographically smaller of kmer and its RC.
    """
    vocab: dict[str, int] = {"[PAD]": 0, "[UNK]": 1, "[MASK]": 2, "[CLS]": 3}
    seen: set[str] = set()
    idx = 4
    for bases in product("ACGT", repeat=k):
        kmer = "".join(bases)
        canonical = min(kmer, _rc_kmer(kmer))
        if canonical not in seen:
            seen.add(canonical)
            vocab[canonical] = idx
            idx += 1
    return vocab


class KmerTokenizer:
    def __init__(self, k: int = 6, stride: int = 3, add_cls: bool = True):
        self.k = k
        self.stride = stride
        self.add_cls = add_cls
        self.vocab = build_canonical_vocab(k)
        self.pad_id = self.vocab["[PAD]"]
        self.unk_id = self.vocab["[UNK]"]
        self.cls_id = self.vocab["[CLS]"]
        self.vocab_size = len(self.vocab)

    def tokenize(self, seq: str) -> list[int]:
        tokens = [self.cls_id] if self.add_cls else []
        for i in range(0, len(seq) - self.k + 1, self.stride):
            kmer = seq[i: i + self.k]
            if "N" in kmer:
                tokens.append(self.unk_id)
            else:
                tokens.append(self.vocab.get(min(kmer, _rc_kmer(kmer)), self.unk_id))
        return tokens

    def __len__(self) -> int:
        return self.vocab_size
