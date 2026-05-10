import torch
import pytest
from src.models.model import _RelativePosBias

def test_relative_pos_bias_shape():
    bias = _RelativePosBias(n_heads=4, max_len=124)
    out = bias(seq_len=124)
    assert out.shape == (4, 124, 124)

def test_relative_pos_bias_shape_short():
    bias = _RelativePosBias(n_heads=4, max_len=124)
    out = bias(seq_len=10)
    assert out.shape == (4, 10, 10)

def test_relative_pos_bias_symmetry():
    bias = _RelativePosBias(n_heads=1, max_len=8)
    with torch.no_grad():
        bias.bias.data = torch.arange(15, dtype=torch.float).unsqueeze(0)
    out = bias(seq_len=8)
    assert out[0, 0, 1].item() == out[0, 1, 2].item()
    assert out[0, 0, 1].item() == out[0, 2, 3].item()
