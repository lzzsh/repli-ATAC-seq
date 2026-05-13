# tests/test_dataset_rt.py
import torch
import pandas as pd
from unittest.mock import patch, MagicMock
from src.data.dataset import RepliSeqDataset, SpeciesConfig

def _make_sp():
    return SpeciesConfig(
        name="mock", fasta="mock.fa", gff3="mock.gff3",
        train_chroms=["chr1"], val_chroms=[], test_chroms=[], species_id=0,
    )

def _make_df(n_bins=600, bin_size=1024):
    classes = ["ES", "MS", "LS", "NR"]
    return pd.DataFrame([{
        "chrom": "chr1",
        "start": i * bin_size,
        "end": (i + 1) * bin_size,
        "RT_class": classes[i % 4],
    } for i in range(n_bins)])

def test_item_rt_labels_shape():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_df(n_bins=400)
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0)
    assert len(ds) > 0
    item = ds[0]
    assert "rt_labels" in item
    assert item["rt_labels"].shape == (896,)
    assert item["rt_labels"].dtype == torch.long

def test_item_no_phase_labels_key():
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 32768 * 4
    mock_genome.fetch.return_value = "A" * 32768
    mock_df = _make_df()
    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=32768, rc_prob=0.0)
    item = ds[0]
    assert "phase_labels" not in item

def test_shift_does_not_produce_negative_win_start():
    """After shift, win_start must be >= 0."""
    sp = _make_sp()
    mock_genome = MagicMock()
    mock_genome.chrom_size.return_value = 196608 * 2
    mock_genome.fetch.return_value = "A" * 196608
    mock_df = _make_df(n_bins=400)

    with patch("src.data.dataset.GenomeSequence", return_value=mock_genome), \
         patch("src.data.dataset.load_labels", return_value=mock_df):
        ds = RepliSeqDataset([sp], "train", window_size=196608, rc_prob=0.0, shift_max=1024)

    # Force win_start=0 so a negative shift would produce a negative coordinate
    ds.samples[0]["win_start"] = 0

    # Patch random.randint to always return -1024 (worst case)
    with patch("src.data.dataset.random.randint", return_value=-1024):
        item = ds[0]

    assert item["one_hot"].shape == (4, 196608)
    assert item["rt_labels"].shape == (896,)
    # Verify fetch was called with start >= 0
    call_args = mock_genome.fetch.call_args
    fetched_start = call_args[0][1]  # positional args: chrom, start, end, window_size
    assert fetched_start >= 0, f"fetch called with negative start: {fetched_start}"
