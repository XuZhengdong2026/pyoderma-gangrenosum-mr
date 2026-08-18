# -*- coding: utf-8 -*-
"""Supplementary Figure S3: cis-pQTL drug-target MR and transcriptome-calibrated
in silico ablation ranking for pyoderma gangrenosum."""
from __future__ import annotations

import csv
import itertools
import os
from math import exp, log

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# MANDATORY font + editable-SVG rules
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42})
mpl.rcParams["font.size"] = 7
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.linewidth"] = 0.8
mpl.rcParams["legend.frameon"] = False

OUT = r"F:\坏疽性脓皮病\outputs\figures"
TABLE = r"F:\坏疽性脓皮病\outputs\tables"

# ----------------------------------------------------------------------------
# Panel A source data (single-instrument Wald ratio; all numbers from analysis)
# target, exposure id, cohort, OR, ci_low, ci_high, p
# ----------------------------------------------------------------------------
mr_rows = [
    ("CXCL8 (IL-8) pQTL\nFolkersen 2020, n = 21,758", "ebi-a-GCST90011994",
     "R5", 1.175, 0.123, 11.22, 0.8884),
    ("CXCL8 (IL-8) pQTL\nFolkersen 2020, n = 21,758", "ebi-a-GCST90011994",
     "R12", 0.870, 0.209, 3.62, 0.8483),
    ("CXCL8 (IL-8) pQTL\nSun 2018, n = 3,301", "prot-a-749",
     "R5", 1.100, 0.227, 5.33, 0.9055),
    ("CXCL8 (IL-8) pQTL\nSun 2018, n = 3,301", "prot-a-749",
     "R12", 0.912, 0.337, 2.47, 0.8569),
    ("IL-1RA (IL1RN) pQTL\nSun 2018, n = 3,301", "prot-a-1504",
     "R5", 0.567, 0.236, 1.36, 0.2053),
    ("IL-1RA (IL1RN) pQTL\nSun 2018, n = 3,301", "prot-a-1504",
     "R12", 0.706, 0.405, 1.23, 0.2177),
    ("TNF pQTL\nSun 2018, n = 3,301", "prot-a-3029",
     "R5", 2.781, 0.603, 12.84, 0.1896),
    ("TNF pQTL\nSun 2018, n = 3,301", "prot-a-3029",
     "R12", 2.618, 1.014, 6.76, 0.0467),
]

# ----------------------------------------------------------------------------
# Panel B source data: hand-set vs transcriptome-calibrated parameter sets
# ----------------------------------------------------------------------------
cal_rows = [
    ("CXCL8 KO", 75.50, 93.12, 48.42, 17.06),
    ("IL1B KO", 100.00, 100.00, 100.00, 100.00),
    ("TNF KO", 28.23, 7.11, 10.70, 0.12),
    ("STAT3 KO", 7.53, 0.10, 2.42, 0.00),
    ("JAK2i", 5.70, 0.06, 1.81, 0.00),
]


def write_source_data() -> None:
    os.makedirs(TABLE, exist_ok=True)
    p1 = os.path.join(TABLE, "Table_S2_pqlt_drug_target_mr.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["target", "exposure_id", "cohort", "OR", "ci_low",
                    "ci_high", "p"])
        w.writerows(mr_rows)
    p2 = os.path.join(TABLE, "Table_S3_calibration_ranking.csv")
    with open(p2, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["intervention", "NEUT_handset_pct", "NEUT_calibrated_pct",
                    "D_handset_pct", "D_calibrated_pct"])
        w.writerows(cal_rows)


def panel_a(ax: plt.Axes) -> None:
    n = len(mr_rows)
    y = np.arange(n)[::-1]
    colors = {"CXCL8": "#0F4D92", "IL-1RA": "#9A4D8E", "TNF": "#B64342"}
    ylabels = []
    for i, (target, eid, cohort, or_, lo, hi, p) in enumerate(mr_rows):
        yi = y[i]
        if target.startswith("CXCL8"):
            c = colors["CXCL8"]
        elif target.startswith("IL-1RA"):
            c = colors["IL-1RA"]
        else:
            c = colors["TNF"]
        # all ORs and CIs are strictly positive, so the log axis is valid
        assert lo > 0 and hi > 0
        lo_c, hi_c = max(lo, 0.05), max(hi, 0.05)
        ax.plot([lo_c, hi_c], [yi, yi], color=c, lw=1.1, zorder=2)
        m = "o" if cohort == "R5" else "s"
        ax.scatter([or_], [yi], marker=m, s=16 if cohort == "R5" else 13,
                   facecolor="white" if cohort == "R12" else c,
                   edgecolor=c, linewidth=0.9, zorder=3)
        txt = f"{or_:.2f} ({lo:.2f}\u2013{hi:.2f})  P = {p:.3f}"
        bold = p < 0.05
        ax.text(hi * 1.05, yi, txt, va="center", ha="left", fontsize=5.9,
                fontweight="bold" if bold else "normal")
        ylabels.append(target if cohort == "R5" else "")
    ax.axvline(1.0, color="#4D4D4D", lw=0.8, ls="--", zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.06, 100)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xticks([0.1, 0.5, 1, 2, 5, 10, 20, 50, 100])
    ax.set_xticklabels(["0.1", "0.5", "1", "2", "5", "10", "20", "50", "100"],
                       fontsize=6)
    ax.set_xlabel("OR per SD increase in genetically proxied protein "
                  "(log scale)", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=6.2)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2.5)
    # legend
    ax.scatter([], [], marker="o", s=16, facecolor="#0F4D92",
               edgecolor="#0F4D92", label="FinnGen R5")
    ax.scatter([], [], marker="s", s=13, facecolor="white",
               edgecolor="#0F4D92", label="FinnGen R12")
    ax.legend(loc="lower right", fontsize=6, frameon=False)


def panel_b(bx: plt.Axes, cx: plt.Axes) -> None:
    labels = [r[0] for r in cal_rows]
    hand_neut = [r[1] for r in cal_rows]
    cal_neut = [r[2] for r in cal_rows]
    hand_d = [r[3] for r in cal_rows]
    cal_d = [r[4] for r in cal_rows]
    x = np.arange(len(labels))
    w = 0.36
    for axx, hand, cal, title in (
            (bx, hand_neut, cal_neut, "Neutrophil infiltration (NEUT)"),
            (cx, hand_d, cal_d, "Ulcer damage (D)")):
        axx.bar(x - w / 2, hand, w, color="#CFCECE", edgecolor="#4D4D4D",
                linewidth=0.5, label="Hand-set parameters")
        axx.bar(x + w / 2, cal, w, color="#42949E", edgecolor="#0F4D92",
                linewidth=0.5, label="Transcriptome-calibrated")
        axx.set_xticks(x)
        axx.set_xticklabels(labels, fontsize=6)
        axx.set_ylim(0, 112)
        axx.set_yticks([0, 25, 50, 75, 100])
        axx.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=6)
        axx.tick_params(axis="x", length=2.5)
        axx.tick_params(axis="y", length=2.5)
        axx.set_ylabel("Reduction vs PG baseline (%)", fontsize=7)
        axx.set_title(title, fontsize=7, pad=2)


def qa_report(fig: plt.Figure, axs) -> None:
    fig.canvas.draw()
    fig_w, fig_h = fig.canvas.get_width_height()
    issues = []
    axes_bbox = [ax.get_window_extent() for ax in axs]
    for ax in axs:
        for t in ax.texts:
            bb = t.get_window_extent()
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fig_w or bb.y1 > fig_h:
                issues.append(f"OUTSIDE-FIGURE: {t.get_text()[:40]!r}")
            for j, other in enumerate(axes_bbox):
                if other is ax.get_window_extent():
                    continue
                if bb.overlaps(other):
                    issues.append(
                        f"CROSS-AXES: {t.get_text()[:30]!r} overlaps axes {j}"
                    )
    for ax in axs:
        texts = [t for t in ax.texts if t.get_text().strip()]
        for a, b in itertools.combinations(texts, 2):
            ba, bb = a.get_window_extent(), b.get_window_extent()
            if not ba.overlaps(bb):
                continue
            inter = ba.intersection(bb)
            area = inter.width * inter.height
            small = min(ba.width * ba.height, bb.width * bb.height)
            if small > 0 and area / small > 0.35:
                issues.append(
                    f"OVERLAP: {a.get_text()[:25]!r} vs {b.get_text()[:25]!r}"
                )
    print("QA issues:", len(issues))
    for x in issues[:50]:
        print(x)


def main() -> None:
    write_source_data()
    fig = plt.figure(figsize=(7.09, 6.1), dpi=600)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.0], wspace=0.28,
                          hspace=0.62, left=0.07, right=0.97, top=0.94,
                          bottom=0.16)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    panel_a(ax_a)
    panel_b(ax_b, ax_c)
    ax_a.text(-0.05, 1.04, "a", transform=ax_a.transAxes, fontsize=8,
              fontweight="bold", va="bottom")
    ax_b.text(-0.04, 1.04, "b", transform=ax_b.transAxes, fontsize=8,
              fontweight="bold", va="bottom")
    ax_b.legend(loc="upper center", fontsize=6, frameon=False, ncol=2,
                bbox_to_anchor=(0.5, -0.30))
    base = os.path.join(OUT, "Figure_S3_pqlt_calibration")
    fig.savefig(base + ".svg", bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(base + ".png", dpi=600, bbox_inches="tight")
    qa_report(fig, [ax_a, ax_b, ax_c])
    plt.close(fig)
    print("Figure_S3_pqlt_calibration written")


if __name__ == "__main__":
    main()
