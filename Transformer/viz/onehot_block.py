"""
One-hot encoded sequence block — matches paper figure style.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch

# ── config ───────────────────────────────────────────────────────────────────
import random
random.seed(7)
# ensure all four bases appear
while True:
    SEQ = ''.join(random.choices('ACGT', k=10))
    if set(SEQ) == {'A','C','G','T'}:
        break
BASE_ORDER = ['A', 'C', 'G', 'T']
BASE_COLOR = {
    'A': '#33A02C',
    'C': '#1F78B4',
    'G': '#FF7F00',
    'T': '#E31A1C',
}

n_cols = len(SEQ)
n_rows = 4
cell  = 0.30     # cell size (data units)
gap   = 0.0      # no gap — cells share borders
pad   = 0.0      # no padding

# grid dimensions
grid_w = n_cols * cell + (n_cols - 1) * gap
grid_h = n_rows * cell + (n_rows - 1) * gap

# figure layout
label_w = 2.2    # left text area width
right_w = 0.5    # right label area
arrow_h = 0.55   # space below grid for arrow
top_h   = 0.25

fig_w = label_w + grid_w + 2 * pad + right_w + 0.3
fig_h = top_h + grid_h + 2 * pad + arrow_h

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis('off')
fig.patch.set_facecolor('white')

# grid origin (bottom-left of first cell)
gx0 = label_w + pad
gy0 = arrow_h + pad

# ── outer border ─────────────────────────────────────────────────────────────
border = patches.Rectangle(
    (gx0, gy0),
    grid_w, grid_h,
    linewidth=2.5, edgecolor='black', facecolor='none', zorder=3
)
ax.add_patch(border)

# ── cells ────────────────────────────────────────────────────────────────────
for row_i, base in enumerate(BASE_ORDER):
    # row 0 = A = top → highest y
    y = gy0 + (n_rows - 1 - row_i) * (cell + gap)
    for col_i, seq_base in enumerate(SEQ):
        x = gx0 + col_i * (cell + gap)

        # colored fill first
        facecolor = BASE_COLOR[seq_base] if seq_base == base else 'white'
        bg = patches.Rectangle(
            (x, y), cell, cell,
            linewidth=0, facecolor=facecolor, zorder=1
        )
        ax.add_patch(bg)

# draw grid lines on top (facecolor=none, just borders)
for row_i in range(n_rows):
    y = gy0 + (n_rows - 1 - row_i) * cell
    for col_i in range(n_cols):
        x = gx0 + col_i * cell
        grid = patches.Rectangle(
            (x, y), cell, cell,
            linewidth=2.0, edgecolor='black', facecolor='none', zorder=3
        )
        ax.add_patch(grid)

# ── right-side base labels ───────────────────────────────────────────────────
for row_i, base in enumerate(BASE_ORDER):
    y = gy0 + (n_rows - 1 - row_i) * (cell + gap) + cell / 2
    x = gx0 + grid_w + pad * 0.8
    ax.text(x, y, base,
            va='center', ha='left',
            fontsize=11, fontweight='bold',
            color=BASE_COLOR[base])

# ── left label ────────────────────────────────────────────────────────────────
mid_y = gy0 + grid_h / 2
ax.text(0.08, mid_y,
        'One-hot encoded\nsequence',
        va='center', ha='left',
        fontsize=11, fontweight='bold',
        color='#5B2D8E',
        linespacing=1.5)

# ── bottom arrow ─────────────────────────────────────────────────────────────
arrow_x = gx0 + grid_w / 2
arrow_y_top = gy0 - pad * 0.6
arrow_y_bot = 0.05

ax.annotate('',
    xy=(arrow_x, arrow_y_bot),
    xytext=(arrow_x, arrow_y_top),
    arrowprops=dict(
        arrowstyle='->', color='#333333',
        lw=1.8,
        mutation_scale=14,
    )
)

plt.tight_layout(pad=0)
for ext in ('png', 'pdf'):
    path = f'/Users/lzz/Documents/GitHub/repli-ATAC-seq/Transformer/viz/onehot_block.{ext}'
    plt.savefig(path, bbox_inches='tight', dpi=300, facecolor='white')
    print('saved:', path)
