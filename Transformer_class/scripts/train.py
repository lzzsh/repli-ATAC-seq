#!/usr/bin/env python3
"""
Single-GPU:
  python scripts/train.py --config src/configs/transformer_wt.yaml

Resume from checkpoint:
  python scripts/train.py --config src/configs/transformer_wt.yaml \
      --resume outputs/basenji2_wt/checkpoints/best_model.pt

Multi-GPU (DDP, e.g. 4 GPUs):
  torchrun --nproc_per_node=4 scripts/train.py --config src/configs/transformer_wt.yaml
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.trainer import train

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--resume", default=None, help="checkpoint path to resume from")
    args = parser.parse_args()
    train(args.config, resume=args.resume)
