import numpy as np
import pandas as pd
import pytest
from src.data.data_utils import load_signals_indexed, SIGNAL_CHANNELS, N_SIGNALS


@pytest.fixture
def sample_df():
    rows = []
    for i in range(50):
        rows.append({
            "chrom": "chr1",
            "start": i * 1024,
            "end": (i + 1) * 1024,
            "G1": float(i),
            "ES": float(i * 2),
            "MS": float(i * 3),
            "LS": float(i * 4),
        })
    return pd.DataFrame(rows)


def test_returns_float32(sample_df):
    query = load_signals_indexed(sample_df)
    result = query("chr1", 0, 28, 1024)
    assert result.dtype == np.float32
    assert result.shape == (28, 4)


def test_signal_values(sample_df):
    query = load_signals_indexed(sample_df)
    result = query("chr1", 0, 1, 1024)
    np.testing.assert_array_almost_equal(result[0], [0.0, 0.0, 0.0, 0.0])
    result2 = query("chr1", 1024, 1, 1024)
    # bin 1: G1=1, ES=2, MS=3, LS=4
    np.testing.assert_array_almost_equal(result2[0], [1.0, 2.0, 3.0, 4.0])


def test_missing_bin_is_nan(sample_df):
    query = load_signals_indexed(sample_df)
    result = query("chr1", 9999 * 1024, 4, 1024)
    assert np.all(np.isnan(result))


def test_center_coordinate_mapping(sample_df):
    query = load_signals_indexed(sample_df)
    result = query("chr1", 0, 1, 1024)
    assert not np.any(np.isnan(result))


def test_mixed_bin_sizes():
    rows = [
        {"chrom": "chr1", "start": 0,    "end": 1024, "G1": 1.0, "ES": 2.0, "MS": 3.0, "LS": 4.0},
        {"chrom": "chr1", "start": 1024, "end": 2048, "G1": 5.0, "ES": 6.0, "MS": 7.0, "LS": 8.0},
        {"chrom": "chr1", "start": 2048, "end": 2560, "G1": 9.0, "ES": 10.0, "MS": 11.0, "LS": 12.0},
    ]
    df = pd.DataFrame(rows)
    query = load_signals_indexed(df)
    result = query("chr1", 2048, 1, 512)
    np.testing.assert_array_almost_equal(result[0], [9.0, 10.0, 11.0, 12.0])


def test_signal_channels_constant():
    assert SIGNAL_CHANNELS == ["G1", "ES", "MS", "LS"]
    assert N_SIGNALS == 4


from src.data.data_utils import one_hot_encode

def test_one_hot_encode_shape():
    arr = one_hot_encode("ACGT")
    assert arr.shape == (4, 4)
    assert arr.dtype == np.float32

def test_one_hot_encode_values():
    arr = one_hot_encode("ACGT")
    np.testing.assert_array_equal(arr[:, 0], [1, 0, 0, 0])  # A
    np.testing.assert_array_equal(arr[:, 1], [0, 1, 0, 0])  # C
    np.testing.assert_array_equal(arr[:, 2], [0, 0, 1, 0])  # G
    np.testing.assert_array_equal(arr[:, 3], [0, 0, 0, 1])  # T

def test_one_hot_encode_unknown_base():
    arr = one_hot_encode("N")
    np.testing.assert_array_equal(arr[:, 0], [0, 0, 0, 0])

def test_one_hot_encode_lowercase():
    arr_lower = one_hot_encode("acgt")
    arr_upper = one_hot_encode("ACGT")
    np.testing.assert_array_equal(arr_lower, arr_upper)

def test_one_hot_encode_long_sequence():
    seq = "ACGT" * 49152  # 196608 bp
    arr = one_hot_encode(seq)
    assert arr.shape == (4, 196608)
    assert arr.sum() == 196608  # exactly one 1 per position
