# tests/test_model_rt.py
import torch
import pytest
from src.models.model import RepliformerModel, RTClassLoss
from src.data.dataset import SpeciesConfig

def _make_sp(name="rice"):
    return SpeciesConfig(name=name, fasta="", gff3="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=0)

def test_model_output_shape():
    model = RepliformerModel([_make_sp()])
    x = torch.zeros(2, 4, 196608)
    out = model(x, head="rice")
    assert "rt_logits" in out
    assert out["rt_logits"].shape == (2, 896, 3)

def test_rt_class_loss():
    criterion = RTClassLoss()
    logits = torch.randn(2, 896, 3)
    labels = torch.randint(0, 3, (2, 896))
    losses = criterion({"rt_logits": logits}, {"rt_labels": labels})
    assert "total" in losses
    assert losses["total"].shape == ()
    assert losses["total"].item() > 0

def test_rt_class_loss_with_weights():
    criterion = RTClassLoss(class_weights=[1.0, 1.0, 1.0])  # 3 classes
    logits = torch.randn(2, 896, 3)
    labels = torch.randint(0, 3, (2, 896))
    losses = criterion({"rt_logits": logits}, {"rt_labels": labels})
    assert losses["total"].shape == ()
