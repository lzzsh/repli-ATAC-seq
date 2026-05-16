import torch
import numpy as np
from unittest.mock import MagicMock, patch
from src.data.dataset import SpeciesConfig


def _make_sp(name, sid):
    return SpeciesConfig(name=name, fasta="", tsv="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=sid)


def _make_batch(n=2):
    return {
        "one_hot":    torch.zeros(n, 4, 196608),
        "rt_signals": torch.rand(n, 896, 4),
        "species_id": torch.zeros(n, dtype=torch.long),
    }


def test_make_loaders_creates_one_loader_per_species():
    cfg = {
        "data": {"input_window_length": 196608},
        "augmentation": {"rc_prob": 0.0, "shift_max": 0},
        "training": {"batch_size": 2},
    }
    mock_ds = MagicMock()
    mock_ds.__len__ = MagicMock(return_value=4)
    with patch("src.trainer.RepliSeqDataset", return_value=mock_ds), \
         patch("src.trainer.DataLoader", return_value=MagicMock()):
        from src.trainer import _make_loaders
        loaders = _make_loaders([_make_sp("rice", 0), _make_sp("maize", 1)],
                                "train", cfg)
    assert set(loaders.keys()) == {"rice", "maize"}


def test_validate_calls_correct_head_per_species():
    called = []

    def fake_forward(one_hot, head):
        called.append(head)
        return {"rt_signals": torch.rand(one_hot.shape[0], 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    criterion = MagicMock(return_value={"total": torch.tensor(0.0), "rt": torch.tensor(0.0)})

    val_loaders = {
        "rice":  [_make_batch(2)],
        "maize": [_make_batch(2)],
    }

    from src.trainer import _validate
    _validate(model, val_loaders, criterion, torch.device("cpu"))
    assert "rice"  in called
    assert "maize" in called


def test_validate_aggregates_metrics_across_species():
    call_count = [0]

    def fake_forward(one_hot, head):
        call_count[0] += 1
        return {"rt_signals": torch.rand(one_hot.shape[0], 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    criterion = MagicMock(return_value={"total": torch.tensor(1.0), "rt": torch.tensor(1.0)})

    val_loaders = {
        "rice":  [_make_batch(2), _make_batch(2)],
        "maize": [_make_batch(2)],
    }

    from src.trainer import _validate
    metrics = _validate(model, val_loaders, criterion, torch.device("cpu"))
    assert call_count[0] == 3
    assert abs(metrics["val_loss_total"] - 1.0) < 1e-5


def test_validate_uses_isnan_mask():
    """_validate must use isnan mask (not != -1) for valid bins."""
    def fake_forward(one_hot, head):
        return {"rt_signals": torch.rand(one_hot.shape[0], 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    criterion = MagicMock(return_value={"total": torch.tensor(0.5), "rt": torch.tensor(0.5)})

    batch = _make_batch(2)
    batch["rt_signals"][0, :10, :] = float("nan")

    val_loaders = {"rice": [batch]}

    from src.trainer import _validate
    metrics = _validate(model, val_loaders, criterion, torch.device("cpu"))
    assert np.isfinite(metrics["val_loss_total"])


def test_validate_reshape_uses_dynamic_n_signals():
    """_validate must not hardcode 4 — reshape must use actual signal dim."""
    n_signals = 2  # deliberately not 4

    def fake_forward(one_hot, head):
        return {"rt_signals": torch.rand(one_hot.shape[0], 896, n_signals)}

    model = MagicMock(side_effect=fake_forward)
    criterion = MagicMock(return_value={"total": torch.tensor(0.0), "rt": torch.tensor(0.0)})

    batch = {
        "one_hot":    torch.zeros(2, 4, 196608),
        "rt_signals": torch.rand(2, 896, n_signals),
        "species_id": torch.zeros(2, dtype=torch.long),
    }
    val_loaders = {"rice": [batch]}

    from src.trainer import _validate
    metrics = _validate(model, val_loaders, criterion, torch.device("cpu"))
    assert "val_loss_total" in metrics
