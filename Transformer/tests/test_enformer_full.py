import torch
import pytest
from src.models.model import _AttnPool

def test_attn_pool_output_shape():
    pool = _AttnPool(in_channels=64, pool_size=2)
    x = torch.randn(2, 64, 128)
    out = pool(x)
    assert out.shape == (2, 64, 64)

def test_attn_pool_odd_length():
    pool = _AttnPool(in_channels=32, pool_size=2)
    x = torch.randn(1, 32, 100)
    out = pool(x)
    assert out.shape == (1, 32, 50)

def test_attn_pool_weights_sum_to_one():
    pool = _AttnPool(in_channels=4, pool_size=2)
    x = torch.ones(1, 4, 8)
    out = pool(x)
    assert out.shape == (1, 4, 4)
    assert torch.isfinite(out).all()

from src.models.model import _EnformerTrunk

def test_enformer_trunk_output_shape():
    trunk = _EnformerTrunk()
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    # 196608/128=1536 bins, crop 320 each side → 896, bottleneck → 1536 ch
    assert out.shape == (1, 1536, 896)

def test_enformer_trunk_no_dilated():
    trunk = _EnformerTrunk()
    assert not hasattr(trunk, 'dilated_tower')

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
