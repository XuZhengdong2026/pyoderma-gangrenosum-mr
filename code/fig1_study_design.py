# -*- coding: utf-8 -*-
"""Figure 1: study design schematic (schematic-led composite, v2.18)."""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 8,
    "axes.linewidth": 0.8,
})

OUT = r"F:\坏疽性脓皮病\outputs\figures"
os.makedirs(OUT, exist_ok=True)

fig, ax = plt.subplots(figsize=(7.5, 5.1))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")


def box(x, y, w, h, text, fc="#F5F5F5", ec="#333333", fs=7.0, bold=False,
        lw=1.1):
    b = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2",
        linewidth=lw, edgecolor=ec, facecolor=fc, mutation_scale=1)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight="bold" if bold else "normal",
            color="#111111")


def arrow(x1, y1, x2, y2, color="#444444", style="-|>", lw=1.2,
          ls="solid"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                        mutation_scale=10, linewidth=lw, color=color,
                        linestyle=ls)
    ax.add_patch(a)


# Row 1: exposures -> instruments -> MR -> PG (R5/R12)
box(1, 84, 15, 12, "GWAS summary statistics\nIBD / CD / UC\n(Liu 2015; de Lange 2017)",
    fc="#EAF1FB")
box(19, 84, 14, 12, "Instruments\nP < 5x10-8\nLD r2 < 0.001 (10 Mb)")
box(36, 84, 14, 12, "Two-sample MR\nIVW primary\nWM / Egger / mode / PRESSO")
box(53, 84, 14, 12, "Steiger + reverse MR\nrelaxed thresholds\nSteiger filtered")
box(70, 84, 17, 12, "Pyoderma gangrenosum\nFinnGen R5 (280 cases)\nR12 replication (703)",
    fc="#FDEBE7", ec="#B71C1C", bold=True)
arrow(16.3, 90, 18.6, 90)
arrow(33.3, 90, 35.6, 90)
arrow(50.3, 90, 52.6, 90)
arrow(67.3, 90, 69.6, 90)

# Mediation branch
box(24, 63, 18, 11, "Two-step mediation\n5 blood-cell traits\nproportions non-significant")
arrow(43, 84, 36, 74.3, color="#B71C1C", ls="dashed")
arrow(33, 69.3, 70, 84, color="#444444", ls="dashed")

# Tissue branch
box(1, 36, 26, 15, "Skin transcriptome\nGSE298908 (58 PG / 45 SS / 10 NN)\n4,285 DEGs + PROGENy\n+ MCP-counter",
    fc="#E8F5E9", ec="#2E7D32")
box(29, 36, 25, 15, "Transcriptome-anchored ODE\nCXCL8 / IL-1 / TNF / STAT3 / JAK2\ngraded + combinatorial ablation",
    fc="#E8F5E9", ec="#2E7D32")
arrow(14, 84, 14, 51.3)
arrow(41, 84, 41, 51.3)

# Drug-target branch
box(57, 36, 20, 15, "cis-pQTL drug-target MR\nCXCL8 / IL-1RA / TNF\nWald-ratio + IVW (R5 / R12)",
    fc="#F3E5F5", ec="#6A1B9A")
box(80, 36, 19, 15, "Colocalization (abf + SuSiE)\nand PheWAS\nnegative / pleiotropy",
    fc="#F3E5F5", ec="#6A1B9A")
arrow(70, 84, 67, 51.3)
arrow(77.5, 43.5, 79.6, 43.5)

# Integration
box(12, 8, 76, 13,
    "Integration (hypothesis-generating): CD-specific causal link + candidate tissue-local\n"
    "CXCL8 / IL-1 programmes; testable treatment hierarchy (CXCR2 / IL-1 / TNF)",
    fc="#FFF8E1", ec="#B26A00", bold=True, fs=7.3)
arrow(20, 36, 32, 21.3)
arrow(50, 36, 50, 21.3)
arrow(76, 36, 64, 21.3)

ax.text(2, 3.2,
        "All genetic analyses were exploratory; nominal P values are reported without global multiple-testing correction.",
        fontsize=6, color="#555555")
ax.text(1, 97.5, "a", fontsize=11, fontweight="bold", color="black")
fig.tight_layout()

base = os.path.join(OUT, "Figure1_study_design")
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{base}.svg", bbox_inches="tight")
fig.savefig(f"{base}.pdf", bbox_inches="tight")
fig.savefig(f"{base}.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("Figure1 saved", fig.get_size_inches())
