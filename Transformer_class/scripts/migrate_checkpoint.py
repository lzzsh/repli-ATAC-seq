# scripts/migrate_checkpoint.py
"""One-time migration: copy old head.{weight,bias} → _heads.<species>.{weight,bias}"""
import argparse
import torch

def migrate(src_path: str, dst_path: str, species_name: str):
    ckpt = torch.load(src_path, map_location="cpu")
    state = ckpt["model"]
    if "head.weight" not in state:
        print("No 'head.weight' found — already migrated or wrong checkpoint.")
        return
    state[f"_heads.{species_name}.weight"] = state.pop("head.weight")
    state[f"_heads.{species_name}.bias"]   = state.pop("head.bias")
    ckpt["model"] = state
    torch.save(ckpt, dst_path)
    print(f"Migrated: head → _heads.{species_name}  →  {dst_path}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("src")
    p.add_argument("dst")
    p.add_argument("--species", required=True)
    args = p.parse_args()
    migrate(args.src, args.dst, args.species)
