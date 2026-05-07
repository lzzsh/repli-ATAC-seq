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
