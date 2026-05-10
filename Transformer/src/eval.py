import torch
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr

from .data.dataset import RepliSeqDataset, load_manifest
from .models.model import Basenji2Model


# ── metrics ───────────────────────────────────────────────────────────────────
def compute_wrt_from_phase_pred(phase_pred: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    # phase_pred is log1p(TPM); invert log1p before TPM-normalizing and computing WRT.
    linear = np.expm1(np.clip(phase_pred, 0, None))  # → TPM scale, non-negative

    def tpm(x):
        return x / (x.sum() + eps) * 1e6

    es = tpm(linear[:, 0])
    ms = tpm(linear[:, 1])
    ls = tpm(linear[:, 2])
    return (0.5 * ms + ls) / (es + ms + ls + eps)


def evaluate_predictions(
    phase_pred: np.ndarray,
    phase_true: np.ndarray,
    wrt_true: np.ndarray,
) -> dict:
    # 支持 [N,4] 和 [N,T,4] 两种shape，统一flatten为 [M,4]
    if phase_pred.ndim == 3:
        N, T, C = phase_pred.shape
        phase_pred = phase_pred.reshape(N * T, C)
        phase_true = phase_true.reshape(N * T, C)
        wrt_true = wrt_true.reshape(N * T)

    valid = np.isfinite(phase_pred).all(axis=1) & np.isfinite(phase_true).all(axis=1)
    phase_pred = phase_pred[valid]
    phase_true = phase_true[valid]
    wrt_true = wrt_true[valid]

    wrt_pred = compute_wrt_from_phase_pred(phase_pred)
    m = {}
    for i, name in enumerate(["ES", "MS", "LS", "G1"]):
        m[f"phase_pearson_{name}"] = pearsonr(phase_pred[:, i], phase_true[:, i])[0]
        m[f"phase_spearman_{name}"] = spearmanr(phase_pred[:, i], phase_true[:, i])[0]
    m["mean_phase_pearson"] = np.mean([m[f"phase_pearson_{n}"] for n in ["ES", "MS", "LS", "G1"]])
    m["mean_phase_spearman"] = np.mean([m[f"phase_spearman_{n}"] for n in ["ES", "MS", "LS", "G1"]])
    m["wrt_pearson"] = pearsonr(wrt_pred, wrt_true)[0]
    m["wrt_spearman"] = spearmanr(wrt_pred, wrt_true)[0]
    return m


# ── shared helpers ────────────────────────────────────────────────────────────
def _load_model(cfg: dict, checkpoint_path: str, device) -> Basenji2Model:
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = Basenji2Model(
        bn_momentum=cfg["model"]["bn_momentum"],
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval()


def _run_inference(model, loader, device):
    pp, pt, wt = [], [], []
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(batch["one_hot"])
            pp.append(out["phase_pred"].cpu().numpy())
            pt_np = batch["phase_labels"].cpu().numpy()
            pt.append(pt_np)
            # derive wrt from labels on the fly
            linear = np.expm1(np.clip(pt_np.reshape(-1, 4), 0, None))
            eps = 1e-6
            wrt = (0.5 * linear[:, 1] + linear[:, 2]) / (linear.sum(axis=1) + eps)
            wt.append(wrt.reshape(pt_np.shape[0], -1))
    return (np.concatenate(pp), np.concatenate(pt), np.concatenate(wt))


# ── WT per-species evaluation ─────────────────────────────────────────────────
def eval_wt(config_path: str, checkpoint_path: str, output_dir: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    species_configs = load_manifest(cfg["data"]["manifest"])
    model = _load_model(cfg, checkpoint_path, device)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in species_configs:
        ds = RepliSeqDataset([sp], "test",
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
    model = _load_model(cfg, checkpoint_path, device)

    target = next(sp for sp in species_configs if sp.name == held_out_species)
    target.test_chroms = target.train_chroms + target.val_chroms + target.test_chroms

    ds = RepliSeqDataset([target], "test",
                         window_size=cfg["data"]["input_window_length"], rc_prob=0.0)
    loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"], shuffle=False, num_workers=4)
    m = evaluate_predictions(*_run_inference(model, loader, device))
    m["held_out_species"] = held_out_species

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([m]).to_csv(out_dir / "cross_species_metrics.tsv", sep="\t", index=False)
    print(m)
