#!/usr/bin/env python3
"""
python scripts/evaluate.py --mode wt \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --output outputs/metrics

python scripts/evaluate.py --mode cross_species \
    --config src/configs/transformer_wt.yaml \
    --checkpoint outputs/basenji2_wt/checkpoints/best_model.pt \
    --output outputs/metrics
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.eval import eval_wt, eval_cross_species

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["wt", "cross_species"], required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--output", default="outputs/metrics")
    args = parser.parse_args()

    if args.mode == "wt":
        eval_wt(args.config, args.checkpoint, args.output)
    else:
        eval_cross_species(args.config, args.checkpoint, args.output)
