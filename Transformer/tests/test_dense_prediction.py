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


from src.models.model import PhaseLoss

def test_phase_loss_dense():
    """PhaseLoss应接受[B,224,4]的pred和true"""
    criterion = PhaseLoss()
    pred = torch.randn(2, 224, 4)
    true = torch.randn(2, 224, 4)
    batch = {"phase_labels": true}
    outputs = {"phase_pred": pred}
    losses = criterion(outputs, batch)
    assert "total" in losses
    assert losses["total"].shape == ()  # scalar
    assert torch.isfinite(losses["total"])


from src.eval import evaluate_predictions, compute_wrt_from_phase_pred

def test_evaluate_predictions_dense():
    """evaluate_predictions应接受[N,224,4]的输入"""
    N = 10
    phase_pred = np.random.rand(N, 224, 4).astype(np.float32)
    phase_true = np.random.rand(N, 224, 4).astype(np.float32)
    wrt_true = np.random.rand(N, 224).astype(np.float32)
    metrics = evaluate_predictions(phase_pred, phase_true, wrt_true)
    assert "phase_pearson_ES" in metrics
    assert "wrt_pearson" in metrics
    assert np.isfinite(metrics["phase_pearson_ES"])


from src.models.model import Basenji2Model, PhaseLoss

def test_trainer_validate_loop_dense():
    """模拟_validate中的数据流，确认dense shape全程兼容"""
    model = Basenji2Model()
    model.eval()
    criterion = PhaseLoss()

    batch = {
        "one_hot": torch.zeros(2, 4, 32768),
        "phase_labels": torch.rand(2, 224, 4),
        "species_id": torch.zeros(2, dtype=torch.long),
    }

    with torch.no_grad():
        out = model(batch["one_hot"])
    losses = criterion(out, batch)
    assert torch.isfinite(losses["total"])

    pp = out["phase_pred"].cpu().numpy()
    pt = batch["phase_labels"].cpu().numpy()
    pt_flat = pt.reshape(-1, 4)
    linear = np.expm1(np.clip(pt_flat, 0, None))
    eps = 1e-6
    wt = (0.5 * linear[:, 1] + linear[:, 2]) / (linear.sum(axis=1) + eps)
    wt = wt.reshape(2, 224)

    metrics = evaluate_predictions(pp, pt, wt)
    assert "phase_pearson_ES" in metrics


def test_end_to_end_forward_and_loss():
    """完整前向传播+loss，确认梯度可以回传"""
    model = Basenji2Model()
    model.train()
    criterion = PhaseLoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)

    batch = {
        "one_hot": torch.zeros(1, 4, 32768),
        "phase_labels": torch.rand(1, 224, 4),
        "species_id": torch.zeros(1, dtype=torch.long),
    }

    out = model(batch["one_hot"])
    losses = criterion(out, batch)
    losses["total"].backward()
    grad = model.trunk.stem.conv.weight.grad
    assert grad is not None
    assert torch.isfinite(grad).all()
    optimizer.step()
