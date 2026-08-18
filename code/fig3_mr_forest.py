# -*- coding: utf-8 -*-
"""Figure 3: MR forest plots (quantitative grid)."""
import math
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

OUT = r"F:\坏疽性脓皮病\outputs\figures"
os.makedirs(OUT, exist_ok=True)

# Panel A: IBD/CD/UC -> PG, methods
panel_a = [
    ("IBD", [
        ("IVW", 0.1505, 0.0855, 0.079),
        ("Weighted median", 0.2183, 0.1311, 0.096),
        ("MR-Egger", -0.3233, 0.2040, 0.115),
    ]),
    ("CD", [
        ("IVW", 0.2132, 0.0826, 0.0099),
        ("Weighted median", 0.2293, 0.1113, 0.039),
        ("MR-Egger", 0.2834, 0.2222, 0.206),
    ]),
    ("UC", [
        ("IVW", -0.0205, 0.0952, 0.830),
        ("Weighted median", -0.0988, 0.1478, 0.504),
        ("MR-Egger", -0.2408, 0.2824, 0.398),
    ]),
]

panel_b = [
    ("Neutrophil count", 0.0722, 0.2704, 0.789),
    ("Neutrophil percentage", 0.4539, 0.2886, 0.116),
    ("Monocyte percentage", 0.0799, 0.2258, 0.724),
    ("White blood cell count", 0.4109, 0.2595, 0.113),
    ("Platelet count", -0.2380, 0.1795, 0.185),
]

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), width_ratios=[1.15, 1.0])
markers = {"IVW": "o", "Weighted median": "s", "MR-Egger": "^"}
colors = {"IVW": "#B71C1C", "Weighted median": "#1565C0", "MR-Egger": "#6A1B9A"}

ax = axes[0]
row = 0
for exposure, methods in panel_a:
    for name, b, se, p in methods:
        lo, hi = b - 1.96 * se, b + 1.96 * se
        ax.errorbar(
            [math.exp(b)], [row], xerr=[[math.exp(b) - math.exp(lo)], [math.exp(hi) - math.exp(b)]],
            fmt=markers[name], color=colors[name], markersize=4.5, capsize=2,
            linewidth=0.9, elinewidth=0.9, zorder=3
        )
        ax.text(3.2, row, f"{math.exp(b):.2f} ({math.exp(lo):.2f}-{math.exp(hi):.2f})", fontsize=5.5, va="center")
        row += 1
    ax.axhline(row - 0.5, color="#CCCCCC", lw=0.5)
ax.set_yticks(range(row))
labels = []
for exposure, methods in panel_a:
    labels += [f"{exposure}: {m[0]}" for m in methods]
ax.set_yticklabels(labels, fontsize=6)
ax.axvline(1, color="black", lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_xlim(0.25, 5.0)
ax.set_xticks([0.25, 0.5, 1, 2, 4])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_xlabel("OR (95% CI) for PG")
ax.set_title("IBD / CD / UC -> PG", fontsize=8)

ax = axes[1]
for i, (name, b, se, p) in enumerate(panel_b):
    lo, hi = b - 1.96 * se, b + 1.96 * se
    ax.errorbar(
        [math.exp(b)], [i], xerr=[[math.exp(b) - math.exp(lo)], [math.exp(hi) - math.exp(b)]],
        fmt="o", color="#1565C0", markersize=4.5, capsize=2,
        linewidth=0.9, elinewidth=0.9, zorder=3
    )
    ax.text(3.0, i, f"{math.exp(b):.2f} ({math.exp(lo):.2f}-{math.exp(hi):.2f})", fontsize=5.5, va="center")
ax.set_yticks(range(len(panel_b)))
ax.set_yticklabels([n for n, *_ in panel_b], fontsize=6)
ax.axvline(1, color="black", lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_xlim(0.25, 5.0)
ax.set_xticks([0.25, 0.5, 1, 2, 4])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_xlabel("OR (95% CI) for PG")
ax.set_title("Blood-cell traits -> PG", fontsize=8)

fig.tight_layout()
base = os.path.join(OUT, "Figure3_MR_forest")
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{base}.svg", bbox_inches="tight")
fig.savefig(f"{base}.pdf", bbox_inches="tight")
fig.savefig(f"{base}.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("Figure3 saved", fig.get_size_inches())
