# tests/test_model_rt.py
import torch
import pytest
from src.models.model import Basenji2Model, RTClassLoss

def test_model_output_shape():
    model = Basenji2Model()
    x = torch.zeros(2, 4, 32768)
    out = model(x)
    assert "rt_logits" in out
    assert out["rt_logits"].shape == (2, 28, 4)

def test_rt_class_loss():
    criterion = RTClassLoss()
    logits = torch.randn(2, 28, 4)
    labels = torch.randint(0, 4, (2, 28))
    losses = criterion({"rt_logits": logits}, {"rt_labels": labels})
    assert "total" in losses
    assert losses["total"].shape == ()
    assert losses["total"].item() > 0
