# tests/test_eval_rt.py
import numpy as np
import pytest
from src.eval import evaluate_predictions

def test_evaluate_predictions_perfect():
    # 3 classes × 3 samples × 28 bins, all predicted correctly
    logits = np.zeros((3, 28, 3))
    labels = np.zeros((3, 28), dtype=np.int64)
    for c in range(3):
        logits[c, :, c] = 10.0
        labels[c, :] = c
    m = evaluate_predictions(logits, labels)
    assert m["macro_f1"] > 0.9
    assert "acc_ES" in m and "acc_MS" in m and "acc_LS" in m
    assert "acc_NR" not in m

def test_evaluate_predictions_flat_input():
    # also accepts [N, 3] logits and [N] labels
    logits = np.zeros((56, 3))
    logits[:, 0] = 10.0
    labels = np.zeros(56, dtype=np.int64)
    m = evaluate_predictions(logits, labels)
    assert m["acc_ES"] == pytest.approx(1.0)
    assert m["macro_f1"] > 0.0

def test_evaluate_predictions_ignore_bins():
    # ignore bins (label=-1) 不参与 metrics
    logits = np.zeros((10, 3))
    logits[:, 0] = 10.0
    labels = np.array([-1, -1, 0, 0, 1, 1, 2, 2, 0, 0], dtype=np.int64)
    m = evaluate_predictions(logits, labels)
    # 只有 label != -1 的 8 个 bin 参与，ES 全对，MS/LS 全错
    assert m["acc_ES"] == pytest.approx(1.0)
    assert m["acc_MS"] == pytest.approx(0.0)
