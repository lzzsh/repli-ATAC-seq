import torch
import pytest
from src.models.model import _EnformerTrunk
from src.models.config_model import RepliformerConfig

def test_enformer_trunk_output_shape():
    trunk = _EnformerTrunk(RepliformerConfig())
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    # 196608/128=1536 tokens, bottleneck → 1536 channels; crop happens in RepliformerModel
    assert out.shape == (1, 1536, 1536)

def test_enformer_trunk_channels():
    trunk = _EnformerTrunk(RepliformerConfig())
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    assert out.shape[1] == 1536

from src.models.model import RepliformerModel, RTSignalLoss
from src.data.dataset import SpeciesConfig

def _make_sp(name="rice"):
    return SpeciesConfig(name=name, fasta="", tsv="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=0)

def test_model_output_shape():
    model = RepliformerModel([_make_sp()])
    x = torch.zeros(1, 4, 196608)
    out = model(x, head="rice")
    assert out["rt_signals"].shape == (1, 896, 4)

def test_model_loss_finite():
    model = RepliformerModel([_make_sp()])
    criterion = RTSignalLoss()
    x = torch.zeros(1, 4, 196608)
    signals = torch.rand(1, 896, 4)
    out = model(x, head="rice")
    loss = criterion(out, {"rt_signals": signals})
    assert torch.isfinite(loss["total"])
    assert loss["total"].item() >= 0

def test_model_param_count():
    model = RepliformerModel([_make_sp()])
    n = sum(p.numel() for p in model.parameters())
    assert 100_000_000 < n < 250_000_000

from src.data.dataset import RepliSeqDataset

def test_dataset_out_bins():
    assert RepliSeqDataset._OUT_BINS == 896

def test_dataset_out_bin_size():
    assert RepliSeqDataset._BIN_SIZE == 128

from unittest.mock import patch as _patch

def test_positional_embed_cached_across_calls():
    """_get_positional_embed must not be called more than once for same T and device."""
    from src.models.model import _EnformerAttention, _get_positional_embed as _orig_fn
    attn = _EnformerAttention(d_model=192, n_heads=4, dim_key=16)
    x = torch.zeros(1, 8, 192)
    call_count = [0]

    def counting_wrapper(*args, **kwargs):
        call_count[0] += 1
        return _orig_fn(*args, **kwargs)

    with _patch("src.models.model._get_positional_embed", side_effect=counting_wrapper):
        attn(x)
        attn(x)

    assert call_count[0] == 1, f"Expected 1 call, got {call_count[0]}"
