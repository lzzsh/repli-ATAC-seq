#!/usr/bin/env python3
"""Replication-timing variant density & purifying selection analysis.

Compares mutation rates and MAF spectra across pure (ES/MS/LS) and mixed
(ESMS/MSLS/ESLS) replicating open chromatin regions.
"""

import sys, time, os
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from intervaltree import IntervalTree
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Paths ──────────────────────────────────────────────────────────────────
PARQUET = "/storage2/liuxiaodongLab/liaozizhuo/Projects/GWAS/rice4k_All_frequency.parquet"
GFF3    = "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3"
OUT_DIR = "/storage2/liuxiaodongLab/liaozizhuo/Projects/GWAS"

CATEGORIES = ["ES", "ESMS", "MS", "MSLS", "LS"]
LABELS     = {"ES": "ES", "MS": "MS", "LS": "LS",
              "ESMS": "ESMS", "MSLS": "MSLS"}
SHORT_LABELS = {"ES": "E", "ESMS": "EM", "MS": "M", "MSLS": "ML", "LS": "L"}
MAF_THRESH = 0.05
MAF_BINS   = [0, 0.01, 0.05, 0.20, 0.51]
MAF_LABELS = ["<1%", "1–5%", "5–20%", "20–50%"]

# ── Nature publication style ────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 7,
    "pdf.fonttype": 42,          # editable TrueType text in PDF
    "svg.fonttype": "none",      # editable <text> in SVG
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.6,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.5,
    "lines.linewidth": 1.0,
})

# Nature colour palette
PALETTE = {
    "blue": "#0F4D92", "teal": "#42949E", "violet": "#9A4D8E",
    "red": "#B64342", "gold": "#E28E2C",
    "neutral_dark": "#4D4D4D", "neutral_mid": "#767676",
    "neutral_light": "#CFCECE", "neutral_pale": "#E8E8E8",
}

# User-specified category colours
CAT_COLORS = {"ES": "#2C5F9E", "ESMS": "#68A0D8", "MS": "#95BE6C",
              "MSLS": "#E4B660", "LS": "#E68364"}
CAT_COLORS_LIST = [CAT_COLORS[c] for c in CATEGORIES]

# MAF bin colours: rare → common (warm to cool)
MAF_COLORS = ["#E8C8C0", "#D4A098", "#A0B8C8", "#688CB8"]

t0 = time.time()

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load RT regions into IntervalTrees
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. Loading RT regions from GFF3 ...")
print("=" * 70)

rt_trees = {cat: defaultdict(IntervalTree) for cat in CATEGORIES}
rt_bp    = {cat: 0 for cat in CATEGORIES}
rt_count = {cat: 0 for cat in CATEGORIES}

with open(GFF3) as f:
    for line in f:
        parts = line.rstrip().split("\t")
        if len(parts) < 4:
            continue
        chrom, start, end, cat = parts[0], int(parts[1]), int(parts[2]), parts[3]
        if cat not in CATEGORIES:
            continue
        if end > start:
            rt_trees[cat][chrom].addi(start, end)
            rt_bp[cat] += end - start
            rt_count[cat] += 1

for cat in CATEGORIES:
    print(f"  {LABELS[cat]:20s}: {rt_count[cat]:>8,} regions, {rt_bp[cat]:>14,} bp, "
          f"avg {rt_bp[cat]/rt_count[cat]:.0f} bp / region")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Merge overlapping intervals & stream parquet to classify variants
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("2. Merging intervals & counting variants per category ...")
print("=" * 70)

# Merge overlapping intervals within each (category, chromosome) to avoid
# double-counting a variant that falls in two overlapping intervals.
rt_merged = {cat: {} for cat in CATEGORIES}
for cat in CATEGORIES:
    for chrom, tree in rt_trees[cat].items():
        if not tree:
            continue
        ivals = sorted(tree, key=lambda x: x.begin)
        merged = []
        for iv in ivals:
            if merged and iv.begin <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], iv.end))
            else:
                merged.append((iv.begin, iv.end))
        rt_merged[cat][chrom] = merged

# Counters: cat -> {all, common, rare, snp, indel, maf_bin: ...}
counts = {cat: {"all": 0, "common": 0, "rare": 0, "snp": 0, "indel": 0,
                **{f"maf_{lb}": 0 for lb in MAF_LABELS}}
          for cat in CATEGORIES}

chroms = [f"chr{i:02d}" for i in range(1, 13)]

for chrom in chroms:
    t_chr = time.time()
    t = pq.read_table(
        PARQUET,
        columns=["position", "MAF", "major_allele", "second_allele"],
        filters=[("chrom", "=", chrom)]
    )
    df = t.to_pandas()
    n = len(df)
    print(f"  {chrom}: {n:>10,} variants", end="", flush=True)

    # Precompute as numpy arrays — vectorized operations below
    pos = df["position"].values
    maf = df["MAF"].values
    is_snp = ((df["major_allele"].str.len() == 1) &
              (df["second_allele"].str.len() == 1)).values
    is_common = maf >= MAF_THRESH
    # digitize: 0→'<1%', 1→'1-5%', 2→'5-20%', 3→'20-50%'
    maf_bin_idx = np.digitize(maf, MAF_BINS) - 1

    for cat in CATEGORIES:
        merged = rt_merged[cat].get(chrom)
        if not merged:
            continue
        # Build boolean mask across all merged intervals
        mask = np.zeros(n, dtype=bool)
        for lo, hi in merged:
            mask |= (pos >= lo) & (pos < hi)

        nc = mask.sum()
        if nc == 0:
            continue

        counts[cat]["all"]    += nc
        counts[cat]["common"] += is_common[mask].sum()
        counts[cat]["rare"]   += nc - is_common[mask].sum()
        counts[cat]["snp"]    += is_snp[mask].sum()
        counts[cat]["indel"]  += nc - is_snp[mask].sum()

        cat_bins = maf_bin_idx[mask]
        for bi, lb in enumerate(MAF_LABELS):
            counts[cat][f"maf_{lb}"] += (cat_bins == bi).sum()

    print(f" ({time.time()-t_chr:.1f}s)")

elapsed = time.time() - t0
print(f"\n  Done streaming in {elapsed:.0f}s")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Variant density summary
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("3. Variant density summary")
print("=" * 70)

rows = []
for cat in CATEGORIES:
    d = counts[cat]
    bp = rt_bp[cat]
    kb = bp / 1000
    rows.append({
        "Category": LABELS[cat],
        "total_bp": bp,
        "n_regions": rt_count[cat],
        "n_variants": d["all"],
        "n_common": d["common"],
        "n_rare": d["rare"],
        "n_snp": d["snp"],
        "n_indel": d["indel"],
        "density_per_kb": d["all"] / kb,
        "density_common_per_kb": d["common"] / kb,
        "density_rare_per_kb": d["rare"] / kb,
        "density_snp_per_kb": d["snp"] / kb,
        "density_indel_per_kb": d["indel"] / kb,
        "prop_common": d["common"] / d["all"] if d["all"] else 0,
        "prop_rare": d["rare"] / d["all"] if d["all"] else 0,
        "prop_snp": d["snp"] / d["all"] if d["all"] else 0,
        "prop_indel": d["indel"] / d["all"] if d["all"] else 0,
    })

density_df = pd.DataFrame(rows)
density_df.to_csv(f"{OUT_DIR}/rt_density_summary.tsv", sep="\t", index=False)
print(density_df.to_string(index=False))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Statistical tests
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("4. Statistical tests")
print("=" * 70)

# 4a. Chi-squared: does variant density differ across ES/MS/LS?
print("\n── 4a. Chi² test: variant count vs bp across ES/MS/LS ──")
observed = np.array([[density_df.iloc[i]["n_variants"],
                      rt_bp[cat] - density_df.iloc[i]["n_variants"]]
                     for i, cat in enumerate(CATEGORIES)])
chi2, p_chi, dof, _ = chi2_contingency(observed)
print(f"  Chi² = {chi2:.1f}, df = {dof}, P = {p_chi:.4e}")
print("  → Densities differ significantly across ES/MS/LS" if p_chi < 0.05
      else "  → No significant difference")

# 4b. Pairwise Fisher exact: common vs rare
print("\n── 4b. Pairwise Fisher exact: common vs rare distribution ──")
for i, c1 in enumerate(CATEGORIES):
    for j, c2 in enumerate(CATEGORIES):
        if i >= j:
            continue
        d1, d2 = counts[c1], counts[c2]
        table = [[d1["common"], d1["rare"]],
                 [d2["common"], d2["rare"]]]
        or_, p = fisher_exact(table, alternative="two-sided")
        print(f"  {c1} vs {c2}: OR={or_:.4f}, P={p:.4e}")

# 4c. Pairwise Fisher exact: SNP vs indel
print("\n── 4c. Pairwise Fisher exact: SNP vs indel distribution ──")
for i, c1 in enumerate(CATEGORIES):
    for j, c2 in enumerate(CATEGORIES):
        if i >= j:
            continue
        d1, d2 = counts[c1], counts[c2]
        table = [[d1["snp"], d1["indel"]],
                 [d2["snp"], d2["indel"]]]
        or_, p = fisher_exact(table, alternative="two-sided")
        print(f"  {c1} vs {c2}: OR={or_:.4f}, P={p:.4e}")

# ══════════════════════════════════════════════════════════════════════════════
# 5. MAF spectrum — purifying selection detection
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("5. MAF spectrum & purifying selection")
print("=" * 70)

maf_rows = []
for cat in CATEGORIES:
    d = counts[cat]
    total = d["all"]
    row = {"Category": LABELS[cat], "total": total}
    for lb in MAF_LABELS:
        v = d[f"maf_{lb}"]
        row[lb] = v
        row[f"{lb}_prop"] = v / total if total else 0
    maf_rows.append(row)

maf_df = pd.DataFrame(maf_rows)
maf_df.to_csv(f"{OUT_DIR}/rt_maf_spectrum.tsv", sep="\t", index=False)

print("\n── 5a. MAF spectrum counts ──")
print(maf_df[["Category", "total"] + MAF_LABELS].to_string(index=False))

print("\n── 5b. MAF bin proportions ──")
for i, cat in enumerate(CATEGORIES):
    props = ", ".join([f"{maf_df.iloc[i][f'{lb}_prop']*100:5.1f}%" for lb in MAF_LABELS])
    print(f"  {LABELS[cat]:20s}: {props}")

# Purifying selection: common(5-50%) vs rare(<5%) ratio
print("\n── 5c. Common/rare ratio (purifying selection index) ──")
all_common = sum(counts[c]["common"] for c in CATEGORIES)
all_rare   = sum(counts[c]["rare"]   for c in CATEGORIES)
bg_cr_ratio = all_common / all_rare  # open-chromatin average
print(f"  Background (ES+MS+LS) common/rare = {bg_cr_ratio:.4f}")

for cat in CATEGORIES:
    d = counts[cat]
    cr = d["common"] / d["rare"] if d["rare"] else 0
    norm = cr / bg_cr_ratio
    # norm < 1 → common depleted relative to bg → purifying selection
    # norm > 1 → common enriched relative to bg → relaxed constraint
    if norm < 0.95:
        signal = "PURIFYING SELECTION (common variants depleted)"
    elif norm > 1.05:
        signal = "RELAXED CONSTRAINT (common variants enriched)"
    else:
        signal = "NEUTRAL"
    print(f"  {LABELS[cat]:20s}: common/rare={cr:.4f}, "
          f"norm={norm:.4f} → {signal}")

# Fisher exact: rare (<1%) vs common (≥5%) enrichment per category vs background
print("\n── 5d. Fisher exact: rare(<1%) vs common(≥5%) vs background ──")
lb_rare = MAF_LABELS[0]
lb_common1 = MAF_LABELS[2]
lb_common2 = MAF_LABELS[3]
bg_rare = sum(counts[c][f"maf_{lb_rare}"] for c in CATEGORIES)
bg_common = sum(counts[c][f"maf_{lb_common1}"] + counts[c][f"maf_{lb_common2}"] for c in CATEGORIES)
print(f"  Background: rare={bg_rare:,}, common={bg_common:,}")

for cat in CATEGORIES:
    d = counts[cat]
    a_rare = d[f"maf_{lb_rare}"]
    a_common = d[f"maf_{lb_common1}"] + d[f"maf_{lb_common2}"]
    b_rare = bg_rare - a_rare
    b_common = bg_common - a_common
    # OR > 1 → category has MORE rare (relative to common) than bg
    # OR < 1 → category has FEWER rare than bg
    table = [[a_rare, a_common], [b_rare, b_common]]
    or_, p = fisher_exact(table, alternative="two-sided")
    direction = "→ more rare than expected" if or_ > 1 else "→ fewer rare than expected"
    print(f"  {LABELS[cat]:20s}: rare={a_rare:,}, common={a_common:,}, "
          f"OR={or_:.4f}, P={p:.4e} {direction}")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Nature-style figure
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("6. Generating Nature-style figure ...")
print("=" * 70)

n_cat = len(CATEGORIES)
short_names = [SHORT_LABELS[c] for c in CATEGORIES]

# ── helpers ─────────────────────────────────────────────────────────────────
def lum(hex_color):
    c = hex_color.lstrip("#")
    return 0.299*int(c[0:2],16) + 0.587*int(c[2:4],16) + 0.114*int(c[4:6],16)

def add_panel_label(ax, label, x=-0.08, y=1.02):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", color="black", ha="left", va="bottom")

# ── figure layout ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
fig.subplots_adjust(wspace=0.42, left=0.08, right=0.96, bottom=0.18, top=0.88)

# ═══ Panel a: Variant density — 3 metric groups × 5 RT categories ══════════
ax = axes[0]
add_panel_label(ax, "a")

all_dens   = np.array([density_df.iloc[i]["density_per_kb"] for i in range(n_cat)])
common_dens = np.array([density_df.iloc[i]["density_common_per_kb"] for i in range(n_cat)])
rare_dens  = np.array([density_df.iloc[i]["density_rare_per_kb"] for i in range(n_cat)])

metric_labels = ["All", "Common", "Rare"]
metric_data   = [all_dens, common_dens, rare_dens]
n_metrics = 3

x = np.arange(n_metrics)
w = 0.15  # bar width per category

for i, cat in enumerate(CATEGORIES):
    vals = [metric_data[m][i] for m in range(n_metrics)]
    offset = (i - (n_cat - 1) / 2) * w
    ax.bar(x + offset, vals, w, color=CAT_COLORS[cat],
           edgecolor="black", lw=0.4, label=SHORT_LABELS[cat])

ax.set_xticks(x)
ax.set_xticklabels(metric_labels)
ax.set_ylabel("Variants per kb")
ax.legend(ncol=5, loc="lower left", bbox_to_anchor=(-0.05, -0.35),
          columnspacing=0.5, handlelength=0.8, handletextpad=0.3)

# chi² annotation
ax.text(0.98, 0.08, f"P = {p_chi:.1e}", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=6, fontstyle="italic", color=PALETTE["neutral_mid"])

# ═══ Panel b: MAF spectrum ══════════════════════════════════════════════════
ax = axes[1]
add_panel_label(ax, "b")

bottom = np.zeros(n_cat)
for i, lb in enumerate(MAF_LABELS):
    vals = np.array([maf_df.iloc[j][f"{lb}_prop"] * 100 for j in range(n_cat)])
    bars = ax.bar(short_names, vals, bottom=bottom, color=MAF_COLORS[i],
                  edgecolor="white", lw=0.3, label=lb)
    # label segments ≥ 8%
    for j, v in enumerate(vals):
        if v >= 8:
            txt_col = "white" if lum(MAF_COLORS[i]) < 128 else "black"
            ax.text(j, bottom[j] + v/2, f"{v:.0f}", ha="center", va="center",
                    fontsize=5.5, fontweight="bold", color=txt_col)
    bottom += vals

ax.set_ylabel("Proportion of variants (%)")
ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(-0.05, -0.35),
          columnspacing=0.6, handlelength=0.7, handletextpad=0.3)

# ═══ Panel c: Common / Rare ratio ═══════════════════════════════════════════
ax = axes[2]
add_panel_label(ax, "c")

cr_vals = np.array([density_df.iloc[i]["n_common"] / density_df.iloc[i]["n_rare"]
                    if density_df.iloc[i]["n_rare"] else 0
                    for i in range(n_cat)])

bars = ax.bar(short_names, cr_vals, color=CAT_COLORS_LIST, edgecolor="white", lw=0.3,
              width=0.65)

# background reference line
ax.axhline(y=bg_cr_ratio, color=PALETTE["neutral_mid"], linestyle="--",
           linewidth=0.8, dashes=(4, 3))
ax.text(0.02, bg_cr_ratio + 0.008, f"avg {bg_cr_ratio:.3f}",
        fontsize=5.5, color=PALETTE["neutral_mid"], va="bottom")

# value + arrow annotations
for bar, val in zip(bars, cr_vals):
    if val < bg_cr_ratio * 0.95:
        tag = " ↓"
    elif val > bg_cr_ratio * 1.05:
        tag = " ↑"
    else:
        tag = ""
    offset = max(cr_vals) * 0.03
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + offset,
            f"{val:.3f}{tag}", ha="center", va="bottom", fontsize=6,
            fontweight="bold", color=PALETTE["neutral_dark"])

ax.set_ylabel("Common / Rare ratio")

# ── export ──────────────────────────────────────────────────────────────────
out_base = f"{OUT_DIR}/rt_variant_analysis"
fig.savefig(f"{out_base}.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{out_base}.svg", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_base}.pdf")
print(f"  Saved: {out_base}.svg")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Done
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"Analysis complete. Total time: {time.time()-t0:.0f}s")
print(f"Outputs:")
print(f"  {OUT_DIR}/rt_density_summary.tsv")
print(f"  {OUT_DIR}/rt_maf_spectrum.tsv")
print(f"  {OUT_DIR}/rt_variant_analysis.pdf")
print(f"  {OUT_DIR}/rt_variant_analysis.svg")
print("=" * 70)
