import torch
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from scipy.stats import pearsonr
from torch.utils.data import DataLoader

from .data.dataset import RepliSeqDataset, load_manifest
from .models.model import RepliformerModel
from .models.config_model import RepliformerConfig

_SIGNAL_CHANNELS = ["G1", "ES", "MS", "LS"]


def evaluate_predictions(
    pred_signals: np.ndarray,
    true_signals: np.ndarray,
) -> dict:
    """
    pred_signals: [N, C] float32
    true_signals: [N, C] float32  (caller has already excluded NaN bins)
    """
    n_channels = pred_signals.shape[1]
    channel_names = _SIGNAL_CHANNELS[:n_channels]
    m = {}
    pearsons = []
    for i, ch in enumerate(channel_names):
        p = pred_signals[:, i]
        t = true_signals[:, i]
        if len(p) > 1 and np.std(t) > 0:
            r, _ = pearsonr(p, t)
        else:
            r = float("nan")
        m[f"pearson_{ch}"] = float(r)
        pearsons.append(r)
    valid_r = [r for r in pearsons if not np.isnan(r)]
    m["mean_pearson"] = float(np.mean(valid_r)) if valid_r else float("nan")
    m["mse"] = float(np.mean((pred_signals - true_signals) ** 2))
    m["mae"] = float(np.mean(np.abs(pred_signals - true_signals)))
    return m


def _load_model(cfg: dict, checkpoint_path: str, device,
                species_configs: list) -> RepliformerModel:
    assert species_configs is not None and len(species_configs) > 0
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = RepliformerModel(
        species_configs=species_configs,
        model_cfg=RepliformerConfig(**cfg.get("model", {})),
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval()


def _run_inference(model, loader, device, head: str):
    pp, pt = [], []
    with torch.no_grad():
        for batch in loader:
            one_hot    = batch["one_hot"].to(device)
            rt_signals = batch["rt_signals"].to(device)
            out = model(one_hot, head=head)
            pp.append(out["rt_signals"].cpu().numpy())
            pt.append(rt_signals.cpu().numpy())
    return np.concatenate(pp), np.concatenate(pt)


def eval_wt(config_path: str, checkpoint_path: str, output_dir: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    species_configs = load_manifest(cfg["data"]["manifest"])
    model = _load_model(cfg, checkpoint_path, device, species_configs)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in species_configs:
        ds = RepliSeqDataset([sp], "test",
                             window_size=cfg["data"]["input_window_length"], rc_prob=0.0)
        loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"],
                            shuffle=False, num_workers=4)
        pred, true = _run_inference(model, loader, device, head=sp.name)
        pred_flat = pred.reshape(-1, pred.shape[-1])
        true_flat = true.reshape(-1, true.shape[-1])
        valid = ~np.isnan(true_flat).any(axis=1)
        m = evaluate_predictions(pred_flat[valid], true_flat[valid])
        m["species"] = sp.name
        rows.append(m)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "wt_test_metrics.tsv", sep="\t", index=False)
    print(df.to_string())


def eval_cross_species(config_path: str, checkpoint_path: str, output_dir: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    species_configs = load_manifest(cfg["data"]["manifest"])
    model = _load_model(cfg, checkpoint_path, device, species_configs)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    loaders = {}
    for sp in species_configs:
        ds = RepliSeqDataset([sp], "test",
                             window_size=cfg["data"]["input_window_length"], rc_prob=0.0)
        loaders[sp.name] = DataLoader(ds, batch_size=cfg["training"]["batch_size"],
                                      shuffle=False, num_workers=4)

    rows = []
    for src_sp in species_configs:
        for tgt_sp in species_configs:
            pred, true = _run_inference(model, loaders[tgt_sp.name], device, head=src_sp.name)
            pred_flat = pred.reshape(-1, pred.shape[-1])
            true_flat = true.reshape(-1, true.shape[-1])
            valid = ~np.isnan(true_flat).any(axis=1)
            m = evaluate_predictions(pred_flat[valid], true_flat[valid])
            m["src_head"]    = src_sp.name
            m["tgt_species"] = tgt_sp.name
            rows.append(m)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "cross_species_metrics.tsv", sep="\t", index=False)
    print(df.to_string())
