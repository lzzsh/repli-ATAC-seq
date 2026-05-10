import torch
import pytest
from src.models.model import _EnformerTrunk

def test_enformer_trunk_output_shape():
    trunk = _EnformerTrunk()
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    # 196608 / 128 = 1536 bins, crop 16 each side → 1504, bottleneck → 3072 channels
    assert out.shape == (1, 3072, 1504)

def test_enformer_trunk_channels():
    trunk = _EnformerTrunk()
    x = torch.zeros(1, 4, 196608)
    out = trunk(x)
    assert out.shape[1] == 3072
