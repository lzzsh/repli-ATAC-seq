import numpy as np
import pandas as pd
import yaml
from unittest.mock import patch, MagicMock
from src.data.dataset import SpeciesConfig


def _make_sp(name, sid):
    return SpeciesConfig(name=name, fasta="f", tsv="g",
                         train_chroms=["c1"], val_chroms=[], test_chroms=["c1"],
                         species_id=sid)


def _mock_inference(model, loader, device, head):
    return np.random.rand(2, 896, 4).astype(np.float32), \
           np.random.rand(2, 896, 4).astype(np.float32)


def _write_cfg(tmp_path):
    cfg = {
        "data": {"manifest": "manifest.yaml", "input_window_length": 196608},
        "model": {"bn_momentum": 0.1},
        "training": {"batch_size": 2},
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.dump(cfg))
    return str(p)


def test_eval_wt_calls_species_head(tmp_path):
    called = []

    def fake_run(model, loader, device, head):
        called.append(head)
        return np.random.rand(1, 896, 4).astype(np.float32), \
               np.random.rand(1, 896, 4).astype(np.float32)

    cfg_path = _write_cfg(tmp_path)
    with patch("src.eval._run_inference", side_effect=fake_run), \
         patch("src.eval._load_model", return_value=MagicMock()), \
         patch("src.eval.load_manifest",
               return_value=[_make_sp("rice", 0), _make_sp("human", 1)]), \
         patch("src.eval.RepliSeqDataset"), \
         patch("src.eval.DataLoader"):
        from src.eval import eval_wt
        eval_wt(cfg_path, "ckpt.pt", str(tmp_path))
    assert called == ["rice", "human"]


def test_eval_cross_species_produces_nxn(tmp_path):
    cfg_path = _write_cfg(tmp_path)
    with patch("src.eval._run_inference", side_effect=_mock_inference), \
         patch("src.eval._load_model", return_value=MagicMock()), \
         patch("src.eval.load_manifest",
               return_value=[_make_sp("rice", 0), _make_sp("human", 1)]), \
         patch("src.eval.RepliSeqDataset"), \
         patch("src.eval.DataLoader"):
        from src.eval import eval_cross_species
        eval_cross_species(cfg_path, "ckpt.pt", str(tmp_path))
    df = pd.read_csv(tmp_path / "cross_species_metrics.tsv", sep="\t")
    assert set(df["src_head"]) == {"rice", "human"}
    assert set(df["tgt_species"]) == {"rice", "human"}
    assert len(df) == 4
