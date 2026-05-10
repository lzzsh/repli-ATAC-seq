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

from src.models.model import _TransformerBlock

def test_transformer_block_output_shape():
    block = _TransformerBlock(d_model=256, n_heads=4, ffn_dim=1024, dropout=0.0)
    x = torch.zeros(2, 124, 256)
    bias = _RelativePosBias(n_heads=4, max_len=124)
    attn_bias = bias(seq_len=124)   # [4, 124, 124]
    out = block(x, attn_bias)
    assert out.shape == (2, 124, 256)

def test_transformer_block_residual():
    block = _TransformerBlock(d_model=4, n_heads=2, ffn_dim=8, dropout=0.0)
    with torch.no_grad():
        for p in block.parameters():
            p.zero_()
    x = torch.ones(1, 3, 4)
    bias = _RelativePosBias(n_heads=2, max_len=3)
    attn_bias = bias(seq_len=3)
    out = block(x, attn_bias)
    assert out.shape == (1, 3, 4)

from src.models.model import Basenji2Model, RTClassLoss

def test_model_output_shape_transformer_head():
    model = Basenji2Model()
    x = torch.zeros(2, 4, 131072)
    out = model(x)
    assert out["rt_logits"].shape == (2, 124, 4)

def test_model_loss_finite():
    model = Basenji2Model()
    criterion = RTClassLoss()
    x = torch.zeros(2, 4, 131072)
    labels = torch.randint(0, 4, (2, 124))
    out = model(x)
    loss = criterion(out, {"rt_labels": labels})
    assert torch.isfinite(loss["total"])
    assert loss["total"].item() > 0
