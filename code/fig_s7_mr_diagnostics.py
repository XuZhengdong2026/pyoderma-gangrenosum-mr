# -*- coding: utf-8 -*-
"""Supplementary Figure S7: CD -> PG MR diagnostics (scatter, funnel,
leave-one-out). Panels rendered by TwoSampleMR; composed here."""
from __future__ import annotations

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42})
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
mpl.rcParams["font.size"] = 7

MR = r"F:\坏疽性脓皮病\outputs\mr"
OUT = r"F:\坏疽性脓皮病\outputs\figures"


def load(path, target_w=900):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    scale = target_w / w
    return im.resize((target_w, int(h * scale)), Image.LANCZOS)


def main():
    scatter = load(os.path.join(MR, "CD_scatter.png"))
    funnel = load(os.path.join(MR, "CD_funnel.png"))
    loo = load(os.path.join(MR, "CD_loo.png"))
    fig = plt.figure(figsize=(7.09, 5.2), dpi=600)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0],
                          height_ratios=[1.0, 1.0], wspace=0.05, hspace=0.08,
                          left=0.02, right=0.98, top=0.97, bottom=0.02)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[:, 1])
    ax_a.imshow(scatter)
    ax_b.imshow(funnel)
    ax_c.imshow(loo)
    for ax in (ax_a, ax_b, ax_c):
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for ax, lab in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
        ax.text(0.015, 0.98, lab, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="left")
    base = os.path.join(OUT, "Figure_S7_mr_diagnostics")
    fig.savefig(base + ".svg", bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(base + ".png", dpi=600, bbox_inches="tight")
    fig.canvas.draw()
    fw, fh = fig.canvas.get_width_height()
    issues = []
    for ax in (ax_a, ax_b, ax_c):
        for t in ax.texts:
            bb = t.get_window_extent()
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fw or bb.y1 > fh:
                issues.append(f"OUTSIDE: {t.get_text()[:30]!r}")
    print("QA issues:", len(issues))
    plt.close(fig)


if __name__ == "__main__":
    main()
