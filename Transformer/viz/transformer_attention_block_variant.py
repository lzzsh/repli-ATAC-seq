"""
Transformer attention block schematic with a distinct visual style.

This variant keeps the same key/query idea but changes the attention paths to
curved flows and uses a blue-green-purple gradient palette.
"""
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path as MplPath
import numpy as np


OUT_DIR = Path(__file__).resolve().parent
OUT_STEM = OUT_DIR / "transformer_attention_block_variant"

FIG_W, FIG_H = 6.8, 3.35
PANEL_EDGE = "#2f3a4a"
TEXT = "#18202f"
MUTED_EDGE = "#bac4d1"
QUERY_EDGE = "#3f2a78"

KEY_GRADIENTS = [
    ("#d9efff", "#7bbcff"),
    ("#bfe7ff", "#3f93ef"),
    ("#94dfd3", "#2fb788"),
    ("#d6f7ba", "#66c97a"),
    ("#c2f0df", "#47bd9d"),
    ("#cdb9ff", "#7b61d1"),
    ("#eadbff", "#ad72f3"),
]
LINE_COLORS = ["#3f93ef", "#2f80ed", "#1ea887", "#56bf76", "#29a587", "#7357ce", "#9c5eea"]
ATTENTION_WEIGHTS = [0.34, 0.46, 0.82, 1.0, 0.56, 0.42, 0.30]


def mix(color_a, color_b, t):
    """Linearly mix two hex colors."""
    a = np.array(mcolors.to_rgb(color_a))
    b = np.array(mcolors.to_rgb(color_b))
    return tuple(a * (1 - t) + b * t)


def gradient_rounded_square(ax, x, y, size, color_left, color_right, edgecolor, lw, radius=0.09, zorder=5):
    """Draw a rounded square filled with a left-to-right gradient."""
    grad = np.linspace(0, 1, 128)
    row = np.array([mix(color_left, color_right, float(t)) for t in grad])
    image = np.tile(row[None, :, :], (128, 1, 1))

    clip = patches.FancyBboxPatch(
        (x, y),
        size,
        size,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=0,
        edgecolor="none",
        facecolor="none",
        transform=ax.transData,
    )
    ax.add_patch(clip)

    im = ax.imshow(
        image,
        extent=(x, x + size, y, y + size),
        origin="lower",
        interpolation="bicubic",
        zorder=zorder,
    )
    im.set_clip_path(clip)

    border = patches.FancyBboxPatch(
        (x, y),
        size,
        size,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor="none",
        joinstyle="round",
        zorder=zorder + 1,
    )
    ax.add_patch(border)
    return border


def rounded_square(ax, x, y, size, facecolor, edgecolor, lw=1.8, radius=0.09, zorder=5):
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


def curved_attention(ax, start, end, color, weight, offset_index, zorder=3):
    """Draw one weighted cubic Bezier attention path."""
    dx = end[0] - start[0]
    side = 1 if dx >= 0 else -1
    curve = 0.22 + 0.035 * abs(offset_index)
    c1 = (start[0] + side * curve, start[1] + 0.34)
    c2 = (end[0] - side * curve * 0.45, end[1] - 0.18)

    path = MplPath(
        [start, c1, c2, end],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    patch = patches.PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        lw=1.1 + 3.3 * weight,
        alpha=0.23 + 0.50 * weight,
        capstyle="round",
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)


def draw_panel():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.set_aspect("equal")
    ax.axis("off")

    panel_x, panel_y = 0.45, 0.35
    panel_w, panel_h = 5.86, 2.56
    panel = patches.FancyBboxPatch(
        (panel_x, panel_y),
        panel_w,
        panel_h,
        boxstyle="round,pad=0.0,rounding_size=0.13",
        linewidth=1.55,
        edgecolor=PANEL_EDGE,
        facecolor="#fbfdff",
        zorder=1,
    )
    ax.add_patch(panel)

    arrow_x = panel_x + panel_w / 2
    ax.annotate(
        "",
        xy=(arrow_x, panel_y + panel_h + 0.02),
        xytext=(arrow_x, FIG_H + 0.12),
        arrowprops=dict(arrowstyle="-|>", lw=1.35, color=PANEL_EDGE, mutation_scale=20),
        clip_on=False,
        zorder=0,
    )
    ax.annotate(
        "",
        xy=(arrow_x, -0.05),
        xytext=(arrow_x, panel_y - 0.02),
        arrowprops=dict(arrowstyle="-|>", lw=1.35, color=PANEL_EDGE, mutation_scale=20),
        clip_on=False,
        zorder=0,
    )

    ax.text(
        panel_x + 0.16,
        panel_y + panel_h - 0.28,
        "Transformer layers (11x)",
        ha="left",
        va="center",
        fontsize=24,
        color=TEXT,
        family="DejaVu Sans",
    )

    n_tokens = 7
    size = 0.47
    gap = 0.115
    start_x = panel_x + 0.55
    key_y = panel_y + 1.54
    query_y = panel_y + 0.29
    xs = [start_x + i * (size + gap) for i in range(n_tokens)]

    ax.text(
        xs[-1] + size + 0.14,
        key_y + size / 2,
        "Key",
        ha="left",
        va="center",
        fontsize=20,
        color=TEXT,
        family="DejaVu Sans",
    )
    ax.text(
        xs[-1] + size + 0.14,
        query_y + size / 2,
        "Query",
        ha="left",
        va="center",
        fontsize=20,
        color=TEXT,
        family="DejaVu Sans",
    )

    query_idx = 3
    query_top_y = query_y + size - 0.01

    for i, x in enumerate(xs):
        rel = i - query_idx
        start = (xs[query_idx] + size / 2 + rel * 0.026, query_top_y)
        end = (x + size / 2, key_y + 0.02)
        curved_attention(ax, start, end, LINE_COLORS[i], ATTENTION_WEIGHTS[i], rel)

    for x, (left, right) in zip(xs, KEY_GRADIENTS):
        gradient_rounded_square(ax, x, key_y, size, left, right, "#1d2430", lw=1.95, radius=0.09, zorder=5)

    for i, x in enumerate(xs):
        if i == query_idx:
            rounded_square(ax, x, query_y, size, "#ffffff", QUERY_EDGE, lw=2.15, radius=0.09, zorder=6)
            inner = patches.Circle(
                (x + size / 2, query_y + size / 2),
                radius=0.065,
                facecolor="#7b61d1",
                edgecolor="none",
                alpha=0.85,
                zorder=7,
            )
            ax.add_patch(inner)
        else:
            rounded_square(ax, x, query_y, size, "#ffffff", MUTED_EDGE, lw=1.6, radius=0.09, zorder=5)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig


if __name__ == "__main__":
    figure = draw_panel()
    for ext in ("png", "pdf", "svg"):
        out_path = OUT_STEM.with_suffix(f".{ext}")
        figure.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.02)
        print(f"saved: {out_path}")
    plt.close(figure)
