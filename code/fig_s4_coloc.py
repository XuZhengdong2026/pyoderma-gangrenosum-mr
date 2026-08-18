# -*- coding: utf-8 -*-
"""Supplementary Figure S4: per-SNP colocalisation posterior probability
(PP.H4) for CD versus plasma cis-pQTL at CXCL8, IL1RN and TNF loci."""
from __future__ import annotations

import csv
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42})
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
mpl.rcParams["font.size"] = 7
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.linewidth"] = 0.8

OUT = r"F:\坏疽性脓皮病\outputs\figures"
INPUT = r"F:\gwas_data\coloc_input"
RESULTS = r"F:\gwas_data\coloc_results"

PANELS = [
    ("a", "CXCL8 (Folkersen 2020)", "coloc_CXCL8_Folkersen2020_vs_CD.csv",
     "pp4_CXCL8_Folkersen2020_vs_CD.csv", 74.55, 74.67),
    ("b", "IL1RN (Sun 2018)", "coloc_IL1RN_Sun2018_vs_CD.csv",
     "pp4_IL1RN_Sun2018_vs_CD.csv", 113.70, 113.95),
    ("c", "IL1RN (Pietzner 2020)", "coloc_IL1RN_Pietzner2020_vs_CD.csv",
     "pp4_IL1RN_Pietzner2020_vs_CD.csv", 113.70, 113.95),
    ("d", "TNF (Sun 2018)", "coloc_TNF_Sun2018_vs_CD.csv",
     "pp4_TNF_Sun2018_vs_CD.csv", 31.10, 31.80),
]


def load(pair_csv, pp_csv):
    pos = {}
    with open(os.path.join(INPUT, pair_csv), encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["snp"]:
                pos[row["snp"]] = float(row["pos"]) / 1e6
    out = []
    with open(os.path.join(RESULTS, pp_csv), encoding="utf-8") as f:
        for row in csv.DictReader(f):
            snp = row["snp"]
            if snp in pos:
                pp = float(row["SNP.PP.H4"])
                out.append((pos[snp], pp, snp))
    out.sort()
    return out


def main():
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.4), dpi=600)
    axes = axes.ravel()
    for ax, (label, title, pair_csv, pp_csv, x0, x1) in zip(axes, PANELS):
        data = load(pair_csv, pp_csv)
        xs = np.array([d[0] for d in data])
        # strictly positive floor (pseudocount) keeps the log axis defined
        ys = np.array([max(d[1], 1e-14) for d in data])
        ax.scatter(xs, ys, s=3, color="#0F4D92", linewidths=0, alpha=0.7)
        imax = int(np.argmax(ys))
        ax.scatter([xs[imax]], [ys[imax]], s=16, facecolor="#B64342",
                   edgecolor="#272727", linewidth=0.6, zorder=3)
        ax.annotate(data[imax][2], (xs[imax], ys[imax]),
                    textcoords="offset points", xytext=(4, 4),
                    fontsize=5.8, color="#272727")
        ax.axvspan(x0, x1, color="#F2F4F7", zorder=0)
        ax.set_yscale("log")
        ax.set_ylim(1e-14, 1e-5)
        ax.set_xlim(min(xs) - 0.02, max(xs) + 0.02)
        ax.set_title(title, fontsize=7, pad=2)
        ax.set_xlabel("Position on chr (Mb)", fontsize=6.5)
        ax.set_ylabel("Per-SNP PP.H4 (log scale)", fontsize=6.5)
        ax.tick_params(axis="x", labelsize=6, length=2.5)
        ax.tick_params(axis="y", labelsize=6, length=2.5)
        ax.text(-0.10, 1.03, label, transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="bottom")
    fig.suptitle("Colocalisation posterior probabilities: CD versus pQTL",
                 fontsize=8, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    base = os.path.join(OUT, "Figure_S4_coloc")
    fig.savefig(base + ".svg", bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(base + ".png", dpi=600, bbox_inches="tight")
    # programmatic QA
    fig.canvas.draw()
    issues = []
    fw, fh = fig.canvas.get_width_height()
    for ax in axes:
        for t in ax.texts:
            bb = t.get_window_extent()
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fw or bb.y1 > fh:
                issues.append(f"OUTSIDE: {t.get_text()[:30]!r}")
    print("QA issues:", len(issues))
    for x in issues[:30]:
        print(x)
    plt.close(fig)


if __name__ == "__main__":
    main()
