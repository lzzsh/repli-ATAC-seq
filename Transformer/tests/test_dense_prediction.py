import numpy as np
import pandas as pd
import pytest
from src.data.data_utils import load_labels_indexed

@pytest.fixture
def sample_df():
    """模拟128bp分辨率的标签DataFrame（已经过load_labels处理）"""
    rows = []
    for i in range(300):
        rows.append({
            "chrom": "chr01",
            "start": i * 128,
            "end": (i + 1) * 128,
            "ES_log1p": float(i) * 0.01,
            "MS_log1p": float(i) * 0.02,
            "LS_log1p": float(i) * 0.03,
            "G1_log1p": float(i) * 0.04,
        })
    return pd.DataFrame(rows)

def test_load_labels_indexed_returns_array(sample_df):
    indexed = load_labels_indexed(sample_df)
    result = indexed("chr01", 0, 224)
    assert result.shape == (224, 4)

def test_load_labels_indexed_values(sample_df):
    indexed = load_labels_indexed(sample_df)
    result = indexed("chr01", 0, 224)
    np.testing.assert_allclose(result[0], [0.0, 0.0, 0.0, 0.0], atol=1e-6)
    np.testing.assert_allclose(result[1, 0], 0.01, atol=1e-6)

def test_load_labels_indexed_missing_bins_are_zero(sample_df):
    """窗口超出标签范围时，缺失bin用0填充"""
    indexed = load_labels_indexed(sample_df)
    result = indexed("chr01", 290, 224)
    assert result.shape == (224, 4)
    assert result[10, 0] == 0.0


from unittest.mock import patch, MagicMock
import pandas as pd
from src.data.dataset import RepliSeqDataset, SpeciesConfig

def _make_mock_species_config():
    return SpeciesConfig(
        name="mock",
        fasta="mock.fa",
        count_tsv="mock.tsv",
        gff3="mock.gff3",
        train_chroms=["chr01"],
        val_chroms=[],
        test_chroms=[],
        species_id=0,
    )

def test_dataset_item_phase_labels_shape():
    """每个样本的phase_labels应为[224, 4]"""
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 32768 * 3
    mock_genome.fetch.return_value = "A" * 32768
    rows = [{"chrom": "chr01", "start": i*128, "end": (i+1)*128,
             "ES_log1p": 0.1, "MS_log1p": 0.2, "LS_log1p": 0.3, "G1_log1p": 0.4,
             "WRT": 0.5, "RT_class": "ES"} for i in range(32768*3//128)]
    mock_df = pd.DataFrame(rows)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=32768, rc_prob=0.0)
    assert len(ds) > 0
    item = ds[0]
    assert item["one_hot"].shape == (4, 32768)
    assert item["phase_labels"].shape == (224, 4)

def test_dataset_no_rt_class_field():
    """改造后不再有rt_class字段"""
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 32768 * 2
    mock_genome.fetch.return_value = "A" * 32768
    rows = [{"chrom": "chr01", "start": i*128, "end": (i+1)*128,
             "ES_log1p": 0.1, "MS_log1p": 0.2, "LS_log1p": 0.3, "G1_log1p": 0.4,
             "WRT": 0.5, "RT_class": "ES"} for i in range(32768*2//128)]
    mock_df = pd.DataFrame(rows)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=32768, rc_prob=0.0)
    item = ds[0]
    assert "rt_class" not in item
    assert "wrt" not in item


import torch
from src.models.model import Basenji2Model

def test_model_output_shape_dense():
    """forward应输出[B, 224, 4]"""
    model = Basenji2Model()
    model.eval()
    x = torch.zeros(2, 4, 32768)
    with torch.no_grad():
        out = model(x)
    assert out["phase_pred"].shape == (2, 224, 4), \
        f"expected (2,224,4), got {out['phase_pred'].shape}"

def test_model_no_center_pooling():
    """不同位置的输出应不同（没有全局pooling）"""
    model = Basenji2Model()
    model.eval()
    x = torch.randn(1, 4, 32768)
    with torch.no_grad():
        out = model(x)
    pred = out["phase_pred"]  # [1, 224, 4]
    assert not torch.allclose(pred[0, 0], pred[0, 112])
