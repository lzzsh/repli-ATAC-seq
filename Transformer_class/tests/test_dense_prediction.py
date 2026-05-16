import numpy as np
import pandas as pd
import torch
import pytest
from unittest.mock import patch, MagicMock
from src.data.data_utils import load_labels_indexed
from src.data.dataset import RepliSeqDataset, SpeciesConfig


@pytest.fixture
def sample_df_1024():
    classes = ["ES", "MS", "LS", "NR"]
    return pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "RT_class": classes[i % 4],
    } for i in range(300)])


def test_load_labels_indexed_returns_array(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr01", 0, 28, 1024)
    assert result.shape == (28,)
    assert result.dtype == np.int64


def test_load_labels_indexed_missing_bins_are_ignore(sample_df_1024):
    query = load_labels_indexed(sample_df_1024)
    result = query("chr01", 9999 * 1024, 28, 1024)
    assert (result == -1).all()  # IGNORE_LABEL=-1


def _make_mock_species_config():
    return SpeciesConfig(
        name="mock", fasta="mock.fa", gff3="mock.gff3",
        train_chroms=["chr01"], val_chroms=[], test_chroms=[], species_id=0,
    )


def test_dataset_item_rt_labels_shape():
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    classes = ["ES", "MS", "LS", "NR"]
    mock_df = pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "RT_class": classes[i % 4],
    } for i in range(196608 * 2 // 1024)])
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    assert len(ds) > 0
    item = ds[0]
    assert item["one_hot"].shape == (4, 196608)
    assert item["rt_labels"].shape == (896,)
    assert item["rt_labels"].dtype == torch.long


def test_dataset_no_phase_labels_key():
    sp = _make_mock_species_config()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 32768 * 2
    mock_genome.fetch.return_value = "A" * 32768
    mock_df = pd.DataFrame([{
        "chrom": "chr01",
        "start": i * 1024,
        "end": (i + 1) * 1024,
        "RT_class": "ES",
    } for i in range(32768 * 2 // 1024)])
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=32768, rc_prob=0.0)
    item = ds[0]
    assert "phase_labels" not in item
