import torch
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import f1_score, roc_auc_score, balanced_accuracy_score

from .data.dataset import RepliSeqDataset, load_manifest
from .tokenization.tokenizer import KmerTokenizer
from .models.model import DNATransformer


# ── metrics ───────────────────────────────────────────────────────────────────
def compute_wrt_from_phase_pred(phase_pred: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    tpm = np.expm1(phase_pred)
    return (0.5 * tpm[:, 1] + tpm[:, 2]) / (tpm.sum(axis=1) + eps)


def evaluate_predictions(
    phase_pred: np.ndarray,
    phase_true: np.ndarray,
    wrt_true: np.ndarray,
    class_pred: np.ndarray,
    class_true: np.ndarray,
) -> dict:
    wrt_pred = compute_wrt_from_phase_pred(phase_pred)
    m = {}
    for i, name in enumerate(["ES", "MS", "LS"]):
        m[f"phase_pearson_{name}"] = pearsonr(phase_pred[:, i], phase_true[:, i])[0]
        m[f"phase_spearman_{name}"] = spearmanr(phase_pred[:, i], phase_true[:, i])[0]
    m["mean_phase_pearson"] = np.mean([m[f"phase_pearson_{n}"] for n in ["ES", "MS", "LS"]])
    m["mean_phase_spearman"] = np.mean([m[f"phase_spearman_{n}"] for n in ["ES", "MS", "LS"]])
    m["wrt_pearson"] = pearsonr(wrt_pred, wrt_true)[0]
    m["wrt_spearman"] = spearmanr(wrt_pred, wrt_true)[0]
    m["macro_f1"] = f1_score(class_true, class_pred, average="macro", zero_division=0)
    m["balanced_accuracy"] = balanced_accuracy_score(class_true, class_pred)
    el = (class_true == 0) | (class_true == 2)
    m["ev_l_auroc"] = float("nan")
    if el.sum() > 10:
        try:
            m["ev_l_auroc"] = roc_auc_score((class_true[el] == 2).astype(int), wrt_pred[el])
        except Exception:
            pass
    return m


# ── shared helpers ────────────────────────────────────────────────────────────
def _load_model(cfg: dict, checkpoint_path: str, n_species: int, vocab_size: int, device) -> DNATransformer:
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = DNATransformer(
        vocab_size=vocab_size, n_species=n_species,
        d_model=cfg["model"]["d_model"], n_layers=cfg["model"]["n_layers"],
        n_heads=cfg["model"]["n_heads"], dim_feedforward=cfg["model"]["dim_feedforward"],
        dropout=0.0, attn_dropout=0.0,
        species_emb_dim=cfg["model"]["species_embedding_dim"],
        head_hidden_dim=cfg["model"]["head_hidden_dim"], head_dropout=0.0,
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval()


def _run_inference(model, loader, device):
    pp, pt, wt, cp, ct = [], [], [], [], []
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(batch["input_ids"], batch["species_id"])
            pp.append(out["phase_pred"].cpu().numpy())
            pt.append(batch["phase_labels"].cpu().numpy())
            wt.append(batch["wrt"].cpu().numpy())
            cp.append(out["class_logits"].argmax(-1).cpu().numpy())
            ct.append(batch["rt_class"].cpu().numpy())
    return (np.concatenate(pp), np.concatenate(pt),
            np.concatenate(wt), np.concatenate(cp), np.concatenate(ct))


# ── WT per-species evaluation ─────────────────────────────────────────────────
def eval_wt(config_path: str, checkpoint_path: str, output_dir: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    species_configs = load_manifest(cfg["data"]["manifest"])
    tokenizer = KmerTokenizer(k=cfg["tokenizer"]["k"], stride=cfg["tokenizer"]["stride"])
    model = _load_model(cfg, checkpoint_path, len(species_configs), tokenizer.vocab_size, device)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in species_configs:
        ds = RepliSeqDataset([sp], "test", tokenizer,
                             window_size=cfg["data"]["input_window_length"], rc_prob=0.0)
        loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"], shuffle=False, num_workers=4)
        m = evaluate_predictions(*_run_inference(model, loader, device))
        m["species"] = sp.name
        rows.append(m)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "wt_test_metrics.tsv", sep="\t", index=False)
    print(df.to_string())


# ── leave-one-species-out evaluation ─────────────────────────────────────────
def eval_cross_species(config_path: str, checkpoint_path: str, held_out_species: str, output_dir: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    species_configs = load_manifest(cfg["data"]["manifest"])
    tokenizer = KmerTokenizer(k=cfg["tokenizer"]["k"], stride=cfg["tokenizer"]["stride"])
    model = _load_model(cfg, checkpoint_path, len(species_configs), tokenizer.vocab_size, device)

    target = next(sp for sp in species_configs if sp.name == held_out_species)
    target.test_chroms = target.train_chroms + target.val_chroms + target.test_chroms

    ds = RepliSeqDataset([target], "test", tokenizer,
                         window_size=cfg["data"]["input_window_length"], rc_prob=0.0)
    loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"], shuffle=False, num_workers=4)
    m = evaluate_predictions(*_run_inference(model, loader, device))
    m["held_out_species"] = held_out_species

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([m]).to_csv(out_dir / "cross_species_metrics.tsv", sep="\t", index=False)
    print(m)
