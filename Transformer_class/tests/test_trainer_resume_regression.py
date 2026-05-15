import ast
from pathlib import Path


def _train_function_tree():
    source = Path("src/trainer.py").read_text()
    module = ast.parse(source)
    for node in module.body:
        if isinstance(node, ast.FunctionDef) and node.name == "train":
            return node
    raise AssertionError("train function not found")


def _first_line(node, predicate):
    lines = [
        child.lineno
        for child in ast.walk(node)
        if hasattr(child, "lineno") and predicate(child)
    ]
    assert lines, "expected node not found"
    return min(lines)


def test_output_dir_is_initialized_before_checkpoint_dir():
    train_fn = _train_function_tree()

    out_dir_assignment = _first_line(
        train_fn,
        lambda node: isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "out_dir" for target in node.targets),
    )
    checkpoint_dir_use = _first_line(
        train_fn,
        lambda node: isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "ckpt_dir" for target in node.targets),
    )

    assert out_dir_assignment < checkpoint_dir_use


def test_gradient_accumulation_is_initialized_before_training_loop():
    train_fn = _train_function_tree()

    accum_assignment = _first_line(
        train_fn,
        lambda node: isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "accum" for target in node.targets),
    )
    accum_use = _first_line(
        train_fn,
        lambda node: isinstance(node, ast.Name)
        and node.id == "accum"
        and isinstance(node.ctx, ast.Load),
    )

    assert accum_assignment < accum_use
