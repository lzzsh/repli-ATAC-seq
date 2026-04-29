import numpy as np
from itertools import product
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import f1_score


# ── feature extractors ────────────────────────────────────────────────────────
def gc_features(seqs: list[str]) -> np.ndarray:
    feats = []
    for s in seqs:
        s = s.upper()
        gc = (s.count("G") + s.count("C")) / max(len(s), 1)
        n_frac = s.count("N") / max(len(s), 1)
        feats.append([gc, n_frac])
    return np.array(feats, dtype=np.float32)


def kmer_features(seqs: list[str], k: int = 4) -> np.ndarray:
    kmers = ["".join(b) for b in product("ACGT", repeat=k)]
    idx = {km: i for i, km in enumerate(kmers)}
    X = np.zeros((len(seqs), len(kmers)), dtype=np.float32)
    for i, seq in enumerate(seqs):
        seq = seq.upper()
        for j in range(len(seq) - k + 1):
            km = seq[j: j + k]
            if km in idx:
                X[i, idx[km]] += 1
        total = X[i].sum()
        if total > 0:
            X[i] /= total
    return X


# ── baseline runners ──────────────────────────────────────────────────────────
def _run_baseline(X_train, X_test, train_wrt, train_class, test_wrt, test_class) -> dict:
    reg = Ridge(alpha=1.0).fit(X_train, train_wrt)
    wrt_pred = reg.predict(X_test)
    clf = LogisticRegression(max_iter=1000, class_weight="balanced", solver="saga").fit(X_train, train_class)
    class_pred = clf.predict(X_test)
    return {
        "wrt_pearson": pearsonr(wrt_pred, test_wrt)[0],
        "wrt_spearman": spearmanr(wrt_pred, test_wrt)[0],
        "macro_f1": f1_score(test_class, class_pred, average="macro", zero_division=0),
    }


def run_gc_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class) -> dict:
    return {"model": "gc", **_run_baseline(
        gc_features(train_seqs), gc_features(test_seqs),
        train_wrt, train_class, test_wrt, test_class,
    )}


def run_kmer_baseline(train_seqs, train_wrt, train_class, test_seqs, test_wrt, test_class, k: int = 4) -> dict:
    return {"model": f"kmer_k{k}", **_run_baseline(
        kmer_features(train_seqs, k), kmer_features(test_seqs, k),
        train_wrt, train_class, test_wrt, test_class,
    )}
