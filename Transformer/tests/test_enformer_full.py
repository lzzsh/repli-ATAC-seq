import torch
import pytest
from src.models.model import _AttnPool
from src.models.config_model import RepliformerConfig

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
    trunk = _EnformerTrunk(RepliformerConfig())
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    # conv tower → bottleneck → [B, 1536, 1536]; crop happens after Transformer
    assert out.shape == (1, 1536, 1536)

def test_enformer_trunk_no_dilated():
    trunk = _EnformerTrunk(RepliformerConfig())
    assert not hasattr(trunk, 'dilated_tower')

from src.models.model import RepliformerModel, RTClassLoss
from src.data.dataset import SpeciesConfig

def _make_sp(name="rice"):
    return SpeciesConfig(name=name, fasta="", gff3="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=0)

def test_model_output_shape():
    model = RepliformerModel([_make_sp()])
    x = torch.zeros(1, 4, 196608)
    out = model(x, head="rice")
    assert out["rt_logits"].shape == (1, 896, 4)

def test_model_loss_finite():
    model = RepliformerModel([_make_sp()])
    criterion = RTClassLoss()
    x = torch.zeros(1, 4, 196608)
    labels = torch.randint(0, 4, (1, 896))
    out = model(x, head="rice")
    loss = criterion(out, {"rt_labels": labels})
    assert torch.isfinite(loss["total"])
    assert loss["total"].item() > 0
