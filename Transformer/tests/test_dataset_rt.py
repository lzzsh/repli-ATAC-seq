# tests/test_dataset_rt.py
import torch
import pandas as pd
from unittest.mock import patch, MagicMock
from src.data.dataset import RepliSeqDataset, SpeciesConfig

def _make_sp():
    return SpeciesConfig(
        name="mock", fasta="mock.fa", count_tsv="mock.tsv", gff3="mock.gff3",
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
    assert item["rt_labels"].shape == (1504,)
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
