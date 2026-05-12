import torch
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from torch.utils.data import DataLoader
from sklearn.metrics import f1_score, confusion_matrix

from .data.dataset import RepliSeqDataset, load_manifest
from .models.model import Basenji2Model

_CLASS_NAMES = ["ES", "MS", "LS", "NR"]


def evaluate_predictions(
    rt_pred_logits: np.ndarray,
    rt_true: np.ndarray,
) -> dict:
    """
    rt_pred_logits: [N, 896, 4] or [N, 4]
    rt_true:        [N, 896] or [N] int64
    """
    if rt_pred_logits.ndim == 3:
        N, T, C = rt_pred_logits.shape
        rt_pred_logits = rt_pred_logits.reshape(N * T, C)
        rt_true = rt_true.reshape(N * T)

    pred_cls = rt_pred_logits.argmax(axis=1)
    m = {}
    for i, name in enumerate(_CLASS_NAMES):
        mask = rt_true == i
        m[f"acc_{name}"] = float((pred_cls[mask] == i).mean()) if mask.any() else float("nan")
    m["macro_f1"] = float(f1_score(rt_true, pred_cls, average="macro", zero_division=0))
    m["overall_acc"] = float((pred_cls == rt_true).mean())
    cm = confusion_matrix(rt_true, pred_cls, labels=list(range(4)))
    for i, name in enumerate(_CLASS_NAMES):
        m[f"cm_row_{name}"] = cm[i].tolist()
    return m


def _load_model(cfg: dict, checkpoint_path: str, device,
                species_configs: list) -> Basenji2Model:
    assert species_configs is not None and len(species_configs) > 0, \
        "species_configs must be provided to _load_model"
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = Basenji2Model(
        species_configs=species_configs,
        bn_momentum=cfg["model"]["bn_momentum"],
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval()


def _run_inference(model, loader, device, head: str):
    pp, pt = [], []
    with torch.no_grad():
        for batch in loader:
            one_hot   = batch["one_hot"].to(device)
            rt_labels = batch["rt_labels"].to(device)
            out = model(one_hot, head=head)
            pp.append(out["rt_logits"].cpu().numpy())
            pt.append(rt_labels.cpu().numpy())
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
        m = evaluate_predictions(*_run_inference(model, loader, device, head=sp.name))
        m["species"] = sp.name
        rows.append(m)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "wt_test_metrics.tsv", sep="\t", index=False)
    print(df.to_string())


def eval_cross_species(config_path: str, checkpoint_path: str, output_dir: str):
    """Evaluate all (src_head × tgt_species) combinations → N×N matrix TSV."""
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
            m = evaluate_predictions(
                *_run_inference(model, loaders[tgt_sp.name], device, head=src_sp.name)
            )
            m["src_head"]    = src_sp.name
            m["tgt_species"] = tgt_sp.name
            rows.append(m)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "cross_species_metrics.tsv", sep="\t", index=False)
    print(df.to_string())
