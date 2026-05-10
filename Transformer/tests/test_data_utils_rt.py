# tests/test_data_utils_rt.py
import numpy as np
import pandas as pd
import pytest
from src.data.data_utils import load_labels_indexed

@pytest.fixture
def sample_df_1024():
    rows = []
    classes = ["ES", "MS", "LS", "NR"]
    for i in range(50):
        rows.append({
            "chrom": "chr1",
            "start": i * 1024,
            "end": (i + 1) * 1024,
            "RT_class": classes[i % 4],
        })
    return pd.DataFrame(rows)

def test_returns_int64(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr1", 0, 28, 1024)
    assert result.dtype == np.int64
    assert result.shape == (28,)

def test_class_values(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr1", 0, 4, 1024)
    # bins 0..3 → ES=0, MS=1, LS=2, NR=3
    np.testing.assert_array_equal(result, [0, 1, 2, 3])

def test_missing_bin_is_nr(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr1", 9999 * 1024, 4, 1024)
    np.testing.assert_array_equal(result, [3, 3, 3, 3])  # NR=3

def test_center_coordinate_mapping(sample_df_1024):
    """Model bin center falls in label bin → correct class."""
    query = load_labels_indexed(sample_df_1024)
    # model bin at genomic pos 0, center = 512, falls in label bin [0,1024) → ES=0
    result = query("chr1", 0, 1, 1024)
    assert result[0] == 0
