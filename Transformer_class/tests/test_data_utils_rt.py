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
    # bins 0..3 → ES=0, MS=1, LS=2, NR=-1 (ignore)
    np.testing.assert_array_equal(result, [0, 1, 2, -1])

def test_missing_bin_is_ignore(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr1", 9999 * 1024, 4, 1024)
    np.testing.assert_array_equal(result, [-1, -1, -1, -1])  # IGNORE_LABEL=-1

def test_center_coordinate_mapping(sample_df_1024):
    """Model bin center falls in label bin → correct class."""
    query = load_labels_indexed(sample_df_1024)
    # model bin at genomic pos 0, center = 512, falls in label bin [0,1024) → ES=0
    result = query("chr1", 0, 1, 1024)
    assert result[0] == 0


def test_mixed_bin_sizes():
    """Last bin is truncated (smaller); its label must still be found."""
    rows = [
        {"chrom": "chr1", "start": 0,    "end": 1024, "RT_class": "ES"},
        {"chrom": "chr1", "start": 1024, "end": 2048, "RT_class": "MS"},
        {"chrom": "chr1", "start": 2048, "end": 2560, "RT_class": "LS"},  # truncated: 512 bp
    ]
    df = pd.DataFrame(rows)
    query = load_labels_indexed(df)
    # model_bin_size=512 → center = 2048 + 256 = 2304, strictly inside [2048, 2560)
    result = query("chr1", 2048, 1, 512)
    assert result[0] == 2  # LS = 2


def test_mixed_bin_sizes_center_outside_uniform_range():
    """Truncated bin whose center does NOT align to the uniform grid."""
    rows = [
        {"chrom": "chr1", "start": 0,    "end": 1024, "RT_class": "ES"},
        {"chrom": "chr1", "start": 1024, "end": 1300, "RT_class": "MS"},  # 276 bp
    ]
    df = pd.DataFrame(rows)
    query = load_labels_indexed(df)
    # model bin center = 1024 + 276//2 = 1162; should map to MS bin [1024, 1300)
    result = query("chr1", 1024, 1, 276)
    assert result[0] == 1  # MS = 1


def test_truncated_bin_center_strictly_inside():
    """Center of model bin falls strictly inside a truncated label bin.

    Old approach: key = (center // label_bin_size) * label_bin_size
    With label_bin_size=1024 (from first bin), center=2200:
      key = (2200 // 1024) * 1024 = 2048 → found in idx_map → would return LS
    But this only works by coincidence. The real failure case is when the
    truncated bin's start does NOT align to the uniform grid.
    """
    rows = [
        {"chrom": "chr1", "start": 0,    "end": 1024, "RT_class": "ES"},
        {"chrom": "chr1", "start": 1024, "end": 2048, "RT_class": "MS"},
        {"chrom": "chr1", "start": 2048, "end": 2700, "RT_class": "LS"},  # 652 bp truncated
    ]
    df = pd.DataFrame(rows)
    query = load_labels_indexed(df)
    # model bin: genomic_start=2048, model_bin_size=652 → center = 2048 + 326 = 2374
    # center 2374 is strictly inside [2048, 2700) → should return LS=2
    result = query("chr1", 2048, 1, 652)
    assert result[0] == 2  # LS = 2
