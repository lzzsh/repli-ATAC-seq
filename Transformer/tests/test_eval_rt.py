# tests/test_eval_rt.py
import numpy as np
import pytest
from src.eval import evaluate_predictions

def test_evaluate_predictions_perfect():
    # 4 samples × 28 bins, all predicted correctly
    logits = np.zeros((4, 28, 4))
    labels = np.zeros((4, 28), dtype=np.int64)
    for c in range(4):
        logits[c % 4, :, c] = 10.0
        labels[c % 4, :] = c
    m = evaluate_predictions(logits, labels)
    assert m["macro_f1"] > 0.9
    assert "acc_ES" in m and "acc_MS" in m and "acc_LS" in m and "acc_NR" in m

def test_evaluate_predictions_flat_input():
    # also accepts [N, 4] logits and [N] labels
    logits = np.zeros((56, 4))
    logits[:, 0] = 10.0
    labels = np.zeros(56, dtype=np.int64)
    m = evaluate_predictions(logits, labels)
    assert m["acc_ES"] == pytest.approx(1.0)
    assert m["macro_f1"] > 0.0
