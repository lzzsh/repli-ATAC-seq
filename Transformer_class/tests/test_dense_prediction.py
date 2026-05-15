import numpy as np
import pandas as pd
import torch
import pytest
from unittest.mock import patch, MagicMock
from src.data.data_utils import load_signals_indexed
from src.data.dataset import RepliSeqDataset, SpeciesConfig


@pytest.fixture
def sample_signal_df():
    return pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "G1": float(i),
        "ES": float(i * 2),
        "MS": float(i * 3),
        "LS": float(i * 4),
    } for i in range(300)])


def test_load_signals_indexed_returns_array(sample_signal_df):
    query = load_signals_indexed(sample_signal_df)
    result = query("chr01", 0, 28, 1024)
    assert result.shape == (28, 4)
    assert result.dtype == np.float32


def test_load_signals_indexed_missing_bins_are_nan(sample_signal_df):
    query = load_signals_indexed(sample_signal_df)
    result = query("chr01", 9999 * 1024, 28, 1024)
    assert np.all(np.isnan(result))


def _make_mock_species_config():
    return SpeciesConfig(
        name="mock", fasta="mock.fa", tsv="mock.tsv",
        train_chroms=["chr01"], val_chroms=[], test_chroms=[], species_id=0,
    )


def test_dataset_item_rt_signals_shape():
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "G1": 1.0, "ES": 2.0, "MS": 3.0, "LS": 4.0,
        "species": "mock",
    } for i in range(196608 * 2 // 1024)])
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    assert len(ds) > 0
    item = ds[0]
    assert item["one_hot"].shape == (4, 196608)
    assert item["rt_signals"].shape == (896, 4)
    assert item["rt_signals"].dtype == torch.float32


def test_dataset_no_rt_labels_key():
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 32768 * 2
    mock_genome.fetch.return_value = "A" * 32768
    mock_df = pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "G1": 1.0, "ES": 2.0, "MS": 3.0, "LS": 4.0,
        "species": "mock",
    } for i in range(32768 * 2 // 1024)])
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=32768, rc_prob=0.0)
    item = ds[0]
    assert "rt_labels" not in item
