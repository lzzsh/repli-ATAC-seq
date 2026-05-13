import torch
import pytest
from src.models.model import _EnformerAttention, _TransformerBlock, _relative_shift

def test_enformer_attention_output_shape():
    attn = _EnformerAttention(d_model=192, n_heads=8, dim_key=64)
    x = torch.zeros(2, 32, 192)
    out = attn(x)
    assert out.shape == (2, 32, 192)

def test_relative_shift_shape():
    # [B, H, T, 2T-1] → [B, H, T, T]
    T = 8
    x = torch.randn(2, 4, T, 2 * T - 1)
    out = _relative_shift(x)
    assert out.shape == (2, 4, T, T)

def test_transformer_block_output_shape():
    # d_model must be divisible by 6 for positional features (num_rel_pos_features = d_model//heads)
    # 192 // 8 = 24, 24 % 6 == 0 ✓
    block = _TransformerBlock(d_model=192, n_heads=8, dim_key=64, dropout=0.0)
    x = torch.zeros(2, 32, 192)
    out = block(x)
    assert out.shape == (2, 32, 192)

def test_transformer_block_residual():
    block = _TransformerBlock(d_model=192, n_heads=8, dim_key=64, dropout=0.0)
    with torch.no_grad():
        for p in block.parameters():
            p.zero_()
    x = torch.ones(1, 16, 192)
    out = block(x)
    assert out.shape == (1, 16, 192)

from src.models.model import RepliformerModel, RTClassLoss
from src.data.dataset import SpeciesConfig

def _make_sp(name="rice"):
    return SpeciesConfig(name=name, fasta="", gff3="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=0)

def test_model_output_shape_transformer_head():
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
