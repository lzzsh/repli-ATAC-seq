"""
Transformer attention block schematic.

Recreates the key/query attention panel from the provided reference image and
saves PNG, PDF, and SVG outputs in this directory.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches


OUT_DIR = Path(__file__).resolve().parent
OUT_STEM = OUT_DIR / "transformer_attention_block"

# Canvas and palette tuned for a paper-figure style.
FIG_W, FIG_H = 6.8, 3.55
BORDER = "#303030"
LIGHT_BORDER = "#a9a9a9"
KEY_EDGE = "#111111"
LINE = "#4fb37a"
KEY_COLORS = [
    "#e3f6df",
    "#bfe8b7",
    "#79cf87",
    "#32a96b",
    "#78cc91",
    "#bce7c4",
    "#e3f6df",
]
QUERY_FILLS = ["#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff"]


def rounded_square(ax, x, y, size, facecolor, edgecolor, lw=2.2, radius=0.095, zorder=3):
    """Draw a rounded square in data coordinates."""
    patch = patches.FancyBboxPatch(
        (x, y),
        size,
        size,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor=facecolor,
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def draw_panel():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.axis("off")

    # Outer rounded transformer-layer container.
    panel_x, panel_y = 0.42, 0.34
    panel_w, panel_h = 5.72, 2.82
    panel = patches.FancyBboxPatch(
        (panel_x, panel_y),
        panel_w,
        panel_h,
        boxstyle="round,pad=0.0,rounding_size=0.12",
        linewidth=1.55,
        edgecolor=BORDER,
        facecolor="white",
        zorder=1,
    )
    ax.add_patch(panel)

    # Incoming and outgoing vertical arrows, clipped by the canvas like the reference.
    arrow_x = panel_x + panel_w / 2
    ax.annotate(
        "",
        xy=(arrow_x, panel_y + panel_h + 0.02),
        xytext=(arrow_x, FIG_H + 0.22),
        arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#222222", mutation_scale=24),
        clip_on=False,
        zorder=0,
    )
    ax.annotate(
        "",
        xy=(arrow_x, -0.08),
        xytext=(arrow_x, panel_y - 0.02),
        arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#222222", mutation_scale=24),
        clip_on=False,
        zorder=0,
    )

    # Title.
    ax.text(
        panel_x + 0.14,
        panel_y + panel_h - 0.28,
        "Transformer layers (11x)",
        ha="left",
        va="center",
        fontsize=25,
        color="#111111",
        family="DejaVu Sans",
    )

    n_tokens = 7
    size = 0.50
    gap = 0.105
    start_x = panel_x + 0.50
    key_y = panel_y + 1.72
    query_y = panel_y + 0.30
    xs = [start_x + i * (size + gap) for i in range(n_tokens)]

    # Label positions.
    ax.text(
        xs[-1] + size + 0.12,
        key_y + size / 2,
        "Key",
        ha="left",
        va="center",
        fontsize=20,
        color="#111111",
        family="DejaVu Sans",
    )
    ax.text(
        xs[-1] + size + 0.12,
        query_y + size / 2,
        "Query",
        ha="left",
        va="center",
        fontsize=20,
        color="#111111",
        family="DejaVu Sans",
    )

    query_idx = 3
    query_center = (xs[query_idx] + size / 2, query_y + size)

    # Attention lines and small key contact dots.
    for i, x in enumerate(xs):
        key_contact = (x + size / 2, key_y + 0.02)
        alpha = 1.0 if i in (2, 3) else 0.48
        lw = 2.1 if i in (2, 3) else 2.0
        ax.plot(
            [query_center[0], key_contact[0]],
            [query_center[1], key_contact[1]],
            color=LINE,
            alpha=alpha,
            lw=lw,
            solid_capstyle="round",
            zorder=2,
        )
        dot = patches.Circle(
            key_contact,
            radius=0.085,
            facecolor=LINE,
            edgecolor=LINE,
            lw=0,
            alpha=alpha,
            zorder=4,
        )
        ax.add_patch(dot)

    # Token boxes above the lines.
    for x, color in zip(xs, KEY_COLORS):
        rounded_square(ax, x, key_y, size, color, KEY_EDGE, lw=2.1, radius=0.095, zorder=5)

    # Query boxes, with the selected query highlighted by a black border.
    for i, (x, fill) in enumerate(zip(xs, QUERY_FILLS)):
        edge = KEY_EDGE if i == query_idx else LIGHT_BORDER
        lw = 2.1 if i == query_idx else 1.9
        rounded_square(ax, x, query_y, size, fill, edge, lw=lw, radius=0.095, zorder=5)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig


if __name__ == "__main__":
    figure = draw_panel()
    for ext in ("png", "pdf", "svg"):
        out_path = OUT_STEM.with_suffix(f".{ext}")
        figure.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.02)
        print(f"saved: {out_path}")
    plt.close(figure)
