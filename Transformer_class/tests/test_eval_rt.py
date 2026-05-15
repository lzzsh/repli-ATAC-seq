import numpy as np
import pytest
from src.eval import evaluate_predictions

_CHANNELS = ["G1", "ES", "MS", "LS"]


def test_evaluate_predictions_perfect_correlation():
    rng = np.random.default_rng(0)
    signals = rng.random((100, 4)).astype(np.float32)
    m = evaluate_predictions(signals, signals)
    for ch in _CHANNELS:
        assert m[f"pearson_{ch}"] == pytest.approx(1.0, abs=1e-5)
    assert m["mean_pearson"] == pytest.approx(1.0, abs=1e-5)


def test_evaluate_predictions_keys():
    rng = np.random.default_rng(1)
    pred = rng.random((50, 4)).astype(np.float32)
    true = rng.random((50, 4)).astype(np.float32)
    m = evaluate_predictions(pred, true)
    for ch in _CHANNELS:
        assert f"pearson_{ch}" in m
    assert "mean_pearson" in m
    assert "mse" in m
    assert "mae" in m


def test_evaluate_predictions_mse_zero_for_perfect():
    rng = np.random.default_rng(2)
    signals = rng.random((50, 4)).astype(np.float32)
    m = evaluate_predictions(signals, signals)
    assert m["mse"] == pytest.approx(0.0, abs=1e-6)
    assert m["mae"] == pytest.approx(0.0, abs=1e-6)


def test_evaluate_predictions_flat_input():
    pred = np.ones((20, 4), dtype=np.float32)
    true = np.ones((20, 4), dtype=np.float32)
    m = evaluate_predictions(pred, true)
    assert "mean_pearson" in m
