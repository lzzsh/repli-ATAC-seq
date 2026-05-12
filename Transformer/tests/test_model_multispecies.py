# tests/test_model_multispecies.py
import torch
import pytest
from src.data.dataset import SpeciesConfig
from src.models.model import Basenji2Model

def _make_configs(names):
    return [SpeciesConfig(name=n, fasta="", gff3="",
                          train_chroms=[], val_chroms=[], test_chroms=[],
                          species_id=i)
            for i, n in enumerate(names)]

def test_heads_created_per_species():
    model = Basenji2Model(_make_configs(["rice", "human"]))
    assert set(model._heads.keys()) == {"rice", "human"}

def test_forward_returns_correct_shape():
    model = Basenji2Model(_make_configs(["rice"]))
    x = torch.zeros(2, 4, 196608)
    out = model(x, head="rice")
    assert out["rt_logits"].shape == (2, 896, 4)

def test_forward_unknown_head_raises():
    model = Basenji2Model(_make_configs(["rice"]))
    x = torch.zeros(1, 4, 196608)
    with pytest.raises(KeyError):
        model(x, head="human")

def test_heads_are_independent():
    model = Basenji2Model(_make_configs(["rice", "human"]))
    rice_params = set(id(p) for p in model._heads["rice"].parameters())
    human_params = set(id(p) for p in model._heads["human"].parameters())
    assert rice_params.isdisjoint(human_params)
