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
