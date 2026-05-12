# tests/test_trainer_routing.py
import torch
from unittest.mock import MagicMock


def _make_batch(species_ids: list[int]):
    B = len(species_ids)
    return {
        "one_hot":    torch.zeros(B, 4, 196608),
        "rt_labels":  torch.zeros(B, 896, dtype=torch.long),
        "species_id": torch.tensor(species_ids, dtype=torch.long),
    }


def test_routing_calls_correct_heads():
    """Each species group must call model with its own head name, exactly once."""
    called_heads = []
    sample_counts = {}

    def fake_forward(one_hot, head):
        called_heads.append(head)
        sample_counts[head] = one_hot.shape[0]
        n = one_hot.shape[0]
        return {"rt_logits": torch.zeros(n, 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    id_to_name = {0: "rice", 1: "human"}

    from src.trainer import _forward_multi_species
    batch = _make_batch([0, 1, 0])  # rice, human, rice
    criterion = MagicMock(return_value={"total": torch.tensor(0.0), "rt": torch.tensor(0.0)})
    _forward_multi_species(model, batch, criterion, id_to_name)
    assert sorted(called_heads) == ["human", "rice"]
    assert called_heads.count("rice") == 1   # called once with 2 samples grouped
    assert called_heads.count("human") == 1  # called once with 1 sample
    # rice was called with 2 samples, human with 1
    assert sample_counts["rice"] == 2
    assert sample_counts["human"] == 1


def test_routing_loss_is_weighted_mean():
    """Loss should be weighted average by group size."""
    def fake_forward(one_hot, head):
        n = one_hot.shape[0]
        return {"rt_logits": torch.zeros(n, 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    id_to_name = {0: "rice", 1: "human"}
    criterion = MagicMock(return_value={"total": torch.tensor(1.0), "rt": torch.tensor(1.0)})

    from src.trainer import _forward_multi_species
    batch = _make_batch([0, 0, 1])  # 2 rice, 1 human
    losses, _ = _forward_multi_species(model, batch, criterion, id_to_name)
    # weighted mean: (2*1.0 + 1*1.0) / 3 = 1.0
    assert abs(losses["total"].item() - 1.0) < 1e-5


def test_routing_logits_shape():
    """Returned logits must be [B, 896, 4] in original order."""
    def fake_forward(one_hot, head):
        n = one_hot.shape[0]
        return {"rt_logits": torch.zeros(n, 896, 4)}

    model = MagicMock(side_effect=fake_forward)
    id_to_name = {0: "rice", 1: "human"}
    criterion = MagicMock(return_value={"total": torch.tensor(0.0), "rt": torch.tensor(0.0)})

    from src.trainer import _forward_multi_species
    batch = _make_batch([0, 1, 0])
    _, logits = _forward_multi_species(model, batch, criterion, id_to_name)
    assert logits.shape == (3, 896, 4)


def test_routing_logits_order():
    """Logits must be placed back in original batch order."""
    def fake_forward(one_hot, head):
        n = one_hot.shape[0]
        val = 1.0 if head == "rice" else 2.0
        return {"rt_logits": torch.full((n, 896, 4), val)}

    model = MagicMock(side_effect=fake_forward)
    id_to_name = {0: "rice", 1: "human"}
    criterion = MagicMock(return_value={"total": torch.tensor(0.0), "rt": torch.tensor(0.0)})

    from src.trainer import _forward_multi_species
    batch = _make_batch([0, 1, 0])  # rice, human, rice
    _, logits = _forward_multi_species(model, batch, criterion, id_to_name)
    assert logits[0, 0, 0].item() == 1.0   # rice
    assert logits[1, 0, 0].item() == 2.0   # human
    assert logits[2, 0, 0].item() == 1.0   # rice
