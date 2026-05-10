import torch
import pytest
from src.models.model import _EnformerTrunk

def test_enformer_trunk_output_shape():
    trunk = _EnformerTrunk()
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    # 196608/128=1536 bins, crop 320 each side → 896, bottleneck → 1536 channels
    assert out.shape == (1, 1536, 896)

def test_enformer_trunk_channels():
    trunk = _EnformerTrunk()
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    assert out.shape[1] == 1536

from src.models.model import Basenji2Model, RTClassLoss

def test_model_output_shape():
    model = Basenji2Model()
    x = torch.zeros(1, 4, 196608)
    out = model(x)
    assert out["rt_logits"].shape == (1, 896, 4)

def test_model_loss_finite():
    model = Basenji2Model()
    criterion = RTClassLoss()
    x = torch.zeros(1, 4, 196608)
    labels = torch.randint(0, 4, (1, 896))
    out = model(x)
    loss = criterion(out, {"rt_labels": labels})
    assert torch.isfinite(loss["total"])
    assert loss["total"].item() > 0

def test_model_param_count():
    model = Basenji2Model()
    n = sum(p.numel() for p in model.parameters())
    assert 100_000_000 < n < 250_000_000

from src.data.dataset import RepliSeqDataset

def test_dataset_out_bins():
    assert RepliSeqDataset._OUT_BINS == 896

def test_dataset_out_bin_size():
    assert RepliSeqDataset._OUT_BIN == 128
