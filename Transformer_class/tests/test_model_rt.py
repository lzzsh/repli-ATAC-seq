import torch
import pytest
from src.models.model import RepliformerModel, RTSignalLoss
from src.data.dataset import SpeciesConfig


def _make_sp(name="rice"):
    return SpeciesConfig(name=name, fasta="", tsv="",
                         train_chroms=[], val_chroms=[], test_chroms=[],
                         species_id=0)


def test_model_output_shape():
    model = RepliformerModel([_make_sp()])
    x = torch.zeros(2, 4, 196608)
    out = model(x, head="rice")
    assert "rt_signals" in out
    assert out["rt_signals"].shape == (2, 896, 4)


def test_model_output_nonnegative():
    """softplus activation ensures predictions >= 0."""
    model = RepliformerModel([_make_sp()])
    x = torch.randn(1, 4, 196608)
    out = model(x, head="rice")
    assert (out["rt_signals"] >= 0).all()


def test_rt_signal_loss_basic():
    criterion = RTSignalLoss()
    pred = torch.rand(2, 896, 4)
    target = torch.rand(2, 896, 4)
    losses = criterion({"rt_signals": pred}, {"rt_signals": target})
    assert "total" in losses and "rt" in losses
    assert losses["total"].shape == ()
    assert losses["total"].item() >= 0


def test_rt_signal_loss_nan_mask():
    """NaN bins in target must be excluded from loss."""
    criterion = RTSignalLoss()
    pred = torch.ones(1, 4, 4)
    target = torch.ones(1, 4, 4)
    target[0, 2, :] = float("nan")
    losses_with_nan = criterion({"rt_signals": pred}, {"rt_signals": target})
    target_no_nan = torch.ones(1, 4, 4)
    losses_no_nan = criterion({"rt_signals": pred}, {"rt_signals": target_no_nan})
    assert torch.isfinite(losses_with_nan["total"])
    assert torch.isfinite(losses_no_nan["total"])


def test_rt_signal_loss_gradient_flows():
    criterion = RTSignalLoss()
    pred = torch.rand(1, 10, 4, requires_grad=True)
    target = torch.rand(1, 10, 4)
    loss = criterion({"rt_signals": pred}, {"rt_signals": target})["total"]
    loss.backward()
    assert pred.grad is not None
