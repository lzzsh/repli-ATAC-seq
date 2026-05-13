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

def test_model_param_count():
    model = RepliformerModel([_make_sp()])
    n = sum(p.numel() for p in model.parameters())
    assert 100_000_000 < n < 250_000_000

from src.data.dataset import RepliSeqDataset

def test_dataset_out_bins():
    assert RepliSeqDataset._OUT_BINS == 896

def test_dataset_out_bin_size():
    assert RepliSeqDataset._OUT_BIN == 128
