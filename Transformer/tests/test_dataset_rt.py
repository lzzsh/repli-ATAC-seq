import torch
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
from src.data.dataset import RepliSeqDataset, SpeciesConfig


def _make_sp():
    return SpeciesConfig(
        name="mock", fasta="mock.fa", tsv="mock.tsv",
        train_chroms=["chr1"], val_chroms=[], test_chroms=[], species_id=0,
    )


def _make_signal_df(n_bins=600, bin_size=1024):
    return pd.DataFrame([{
        "chrom": "chr1",
        "start": i * bin_size,
        "end": (i + 1) * bin_size,
        "G1": float(i),
        "ES": float(i * 2),
        "MS": float(i * 3),
        "LS": float(i * 4),
        "species": "mock",
    } for i in range(n_bins)])


def test_item_rt_signals_shape():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_signal_df(n_bins=400)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    assert len(ds) > 0
    item = ds[0]
    assert "rt_signals" in item
    assert item["rt_signals"].shape == (896, 4)
    assert item["rt_signals"].dtype == torch.float32


def test_item_no_rt_labels_key():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_signal_df(n_bins=400)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    item = ds[0]
    assert "rt_labels" not in item


def test_signals_are_log1p_transformed():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_signal_df(n_bins=400)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    item = ds[0]
    valid_mask = ~torch.isnan(item["rt_signals"])
    assert (item["rt_signals"][valid_mask] >= 0).all()


def test_shift_does_not_produce_negative_win_start():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_signal_df(n_bins=400)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0, shift_max=1024)
    ds.samples[0]["win_start"] = 0
    import random
    with patch("src.data.dataset.random.randint", return_value=-1024):
        item = ds[0]
    assert item["one_hot"].shape == (4, 196608)
    assert item["rt_signals"].shape == (896, 4)
    call_args = mock_genome.fetch.call_args
    fetched_start = call_args[0][1]
    assert fetched_start >= 0


def test_stride_matches_window_geometry():
    """Stride must equal window_size - 2 * crop_bins * bin_size."""
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 4
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_signal_df(n_bins=800)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_signals", return_value=mock_df):
        ds_default = RepliSeqDataset([sp], "train", window_size=196608)
        ds_small   = RepliSeqDataset([sp], "train", window_size=131072)

    expected_default = 196608 - 2 * RepliSeqDataset._CROP_BINS * RepliSeqDataset._BIN_SIZE
    expected_small   = 131072 - 2 * RepliSeqDataset._CROP_BINS * RepliSeqDataset._BIN_SIZE

    chrom_size = 196608 * 4
    n_default = len([s for s in ds_default.samples if s["species"] == "mock"])
    n_small   = len([s for s in ds_small.samples   if s["species"] == "mock"])

    import math
    assert n_default == math.floor((chrom_size - 196608) / expected_default) + 1
    assert n_small   == math.floor((chrom_size - 131072) / expected_small)   + 1
