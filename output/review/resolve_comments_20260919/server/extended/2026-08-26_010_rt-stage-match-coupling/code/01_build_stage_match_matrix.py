#!/usr/bin/env python3
"""Build the WT RT-class by ATAC-stage long matrix from task 009 OCR metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


STAGES = ("ES", "MS", "LS")
ALLELES = ("cpp8_1", "cpp8_3")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ocr-metrics", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    output = Path(args.output).resolve()
    task_root = Path(__file__).resolve().parents[1]
    if task_root not in output.parents:
        raise ValueError(f"Refusing output outside task root: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    ocr = pd.read_csv(args.ocr_metrics, sep="\t")
    phase_cols = ["WT_ES_norm", "WT_MS_norm", "WT_LS_norm"]
    required = {"chr", "start", "end", "block", "length", "WT_wrt", *phase_cols}
    for allele in ALLELES:
        required.add(f"delta_rt_{allele}")
        for stage in STAGES:
            required.update({f"delta_atac_{allele}_{stage.lower()}", f"wt_{stage.lower()}_log_atac"})
    missing = required - set(ocr.columns)
    if missing:
        raise ValueError(f"Missing task-009 columns: {sorted(missing)}")

    signals = ocr[phase_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    total = signals.sum(axis=1)
    proportions = np.divide(signals, total[:, None], out=np.full_like(signals, np.nan), where=total[:, None] > 0)
    order = np.sort(proportions, axis=1)
    confidence = order[:, -1] - order[:, -2]
    class_index = np.argmax(np.where(np.isfinite(proportions), proportions, -np.inf), axis=1)
    class_names = np.array(["E", "M", "L"])[class_index]
    valid_phase = np.isfinite(proportions).all(axis=1) & (total > 0)

    base = ocr[["chr", "start", "end", "block", "length", "WT_wrt"]].copy()
    base["wt_rt_class"] = class_names
    base["rt_class_index"] = class_index
    base["rt_confidence"] = confidence
    for i, stage in enumerate(STAGES):
        base[f"wt_{stage.lower()}_phase_proportion"] = proportions[:, i]
    base = base[valid_phase].copy()

    frames = []
    for allele in ALLELES:
        for stage_index, stage in enumerate(STAGES):
            sl = stage.lower()
            d = base.copy()
            d["allele"] = allele
            d["atac_stage"] = stage
            d["atac_stage_index"] = stage_index
            d["lag"] = stage_index - d["rt_class_index"]
            d["cell"] = d["wt_rt_class"] + "_" + stage
            d["delta_rt"] = pd.to_numeric(ocr.loc[d.index, f"delta_rt_{allele}"], errors="coerce")
            d["delta_atac"] = pd.to_numeric(ocr.loc[d.index, f"delta_atac_{allele}_{sl}"], errors="coerce")
            d["wt_log_atac"] = pd.to_numeric(ocr.loc[d.index, f"wt_{sl}_log_atac"], errors="coerce")
            frames.append(d)
    long = pd.concat(frames, ignore_index=True)
    long.to_csv(output, sep="\t", index=False)
    print("OCRs with WT phase class:", base.shape[0])
    print("Long observations:", long.shape[0])
    print(base["wt_rt_class"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
