#!/usr/bin/env python3
"""
python scripts/virtual_knockout.py saturation --config ... --checkpoint ... \
    --fasta ... --species rice --species_id 0 --n_species 3 \
    --chrom chr1 --start 1000000 --end 1008192 \
    --output outputs/interpretation/saturation.tsv

python scripts/virtual_knockout.py motif ... --motif_start 4000 --motif_end 4020 --mut_mode scramble
python scripts/virtual_knockout.py flank ... --motif_start 4000 --motif_end 4020
python scripts/virtual_knockout.py ig    ... --output outputs/interpretation/ig.npy
"""
import argparse
import sys
import yaml
import torch
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data.data_utils import GenomeSequence
from src.tokenization.tokenizer import KmerTokenizer
from src.models.model import DNATransformer
from src.interpret import (
    saturation_mutagenesis, importance_per_position,
    motif_mutation, flank_control, integrated_gradients,
)


def _load(cfg, ckpt_path, n_species, device):
    tokenizer = KmerTokenizer(k=cfg["tokenizer"]["k"], stride=cfg["tokenizer"]["stride"])
    ckpt = torch.load(ckpt_path, map_location=device)
    model = DNATransformer(
        vocab_size=tokenizer.vocab_size, n_species=n_species,
        d_model=cfg["model"]["d_model"], n_layers=cfg["model"]["n_layers"],
        n_heads=cfg["model"]["n_heads"], dim_feedforward=cfg["model"]["dim_feedforward"],
        dropout=0.0, attn_dropout=0.0,
        species_emb_dim=cfg["model"]["species_embedding_dim"],
        head_hidden_dim=cfg["model"]["head_hidden_dim"], head_dropout=0.0,
    )
    model.load_state_dict(ckpt["model"])
    return model.to(device).eval(), tokenizer


def _predict_fn(model, tokenizer, species_id_val, device):
    def fn(seq: str) -> np.ndarray:
        ids = torch.tensor([tokenizer.tokenize(seq)], dtype=torch.long).to(device)
        sp = torch.tensor([species_id_val], dtype=torch.long).to(device)
        with torch.no_grad():
            return model(ids, sp)["phase_pred"].squeeze(0).cpu().numpy()
    return fn


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["saturation", "motif", "flank", "ig"])
    p.add_argument("--config", required=True)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--species_id", type=int, required=True)
    p.add_argument("--n_species", type=int, default=3)
    p.add_argument("--chrom", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--motif_start", type=int, default=None)
    p.add_argument("--motif_end", type=int, default=None)
    p.add_argument("--mut_mode", default="scramble", choices=["scramble", "gc_matched_random"])
    p.add_argument("--output", required=True)
    args = p.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, tokenizer = _load(cfg, args.checkpoint, args.n_species, device)
    seq = GenomeSequence(args.fasta).fetch(args.chrom, args.start, args.end,
                                           cfg["data"]["input_window_length"])
    predict_fn = _predict_fn(model, tokenizer, args.species_id, device)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "saturation":
        results = saturation_mutagenesis(seq, predict_fn)
        pd.DataFrame(results).to_csv(out, sep="\t", index=False)
        np.save(str(out).replace(".tsv", "_importance.npy"), importance_per_position(results))
        print(f"Saved {len(results)} mutations → {out}")

    elif args.mode == "motif":
        assert args.motif_start is not None
        r = motif_mutation(seq, args.motif_start, args.motif_end, predict_fn, args.mut_mode)
        pd.DataFrame([r]).to_csv(out, sep="\t", index=False)
        print(r)

    elif args.mode == "flank":
        assert args.motif_start is not None
        results = flank_control(seq, args.motif_start, args.motif_end, predict_fn)
        pd.DataFrame(results).to_csv(out, sep="\t", index=False)
        print(pd.DataFrame(results).to_string())

    elif args.mode == "ig":
        ids = torch.tensor([tokenizer.tokenize(seq)], dtype=torch.long)
        sp = torch.tensor([args.species_id], dtype=torch.long)
        attrs = integrated_gradients(model, ids, sp, target_idx=3, device=str(device))
        np.save(str(out).replace(".tsv", ".npy"), attrs)
        print(f"IG shape: {attrs.shape} → {out}")


if __name__ == "__main__":
    main()
