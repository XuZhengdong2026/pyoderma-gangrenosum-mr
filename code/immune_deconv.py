# -*- coding: utf-8 -*-
"""MCP-counter immune/stromal cell abundance estimates on GSE298908
(log2 median-of-ratios normalised counts; marker means per population)."""
from __future__ import annotations

import csv
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42})
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
mpl.rcParams["font.size"] = 7
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.linewidth"] = 0.8

COUNTS = r"F:\tx_data\GSE298908_counts.txt"
MARKERS = r"F:\gwas_data\mcpcounter\MCPcounter-master\Signatures\genes.txt"
OUT = r"F:\坏疽性脓皮病\outputs\figures"
TAB = r"F:\坏疽性脓皮病\outputs\tables"


def read_counts(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        samples = header
        genes, rows = [], []
        for line in reader:
            if not line or not line[0]:
                continue
            genes.append(line[0])
            rows.append([float(x) for x in line[1:1 + len(samples)]])
    return np.array(rows), genes, samples


def group_of(col):
    if "Pyoderma_Gangrenosum" in col:
        return "PG"
    if "Sweet" in col:
        return "SS"
    if re.search(r"\.NN\.", col):
        return "NN"
    return "NA"


def deseq_log2(counts):
    gm = np.exp(np.mean(np.log(counts + 0.5), axis=1))
    gm[gm <= 0] = np.nan
    ratio = counts / gm[:, None]
    with np.errstate(invalid="ignore"):
        sf = np.nanmedian(np.where(counts > 0, ratio, np.nan), axis=0)
    sf = np.where(np.isfinite(sf) & (sf > 0), sf, np.nanmedian(sf))
    # +1 pseudocount keeps the log transform defined for zero counts
    return np.log2(counts / sf + 1)


def load_markers(path):
    df = pd.read_csv(path, sep="\t")
    out = {}
    for pop, g in df.groupby("Cell population"):
        out[pop] = set(g["HUGO symbols"].str.upper())
    return out


def main():
    counts, genes, samples = read_counts(COUNTS)
    expr = deseq_log2(counts)
    groups = np.array([group_of(s) for s in samples])
    markers = load_markers(MARKERS)
    gset = set(genes)
    scores = {}
    used = {}
    for pop, mk in markers.items():
        m = sorted(mk & gset)
        idx = [genes.index(g) for g in m]
        scores[pop] = expr[idx, :].mean(axis=0)
        used[pop] = len(m)
    score_df = pd.DataFrame(scores, index=samples)
    score_df.insert(0, "group", groups)
    score_df.to_csv(os.path.join(TAB, "Table_S6_mcpcounter_scores.csv"),
                    encoding="utf-8")

    rows = []
    for pop in score_df.columns[1:]:
        for ref, case in (("NN", "PG"), ("SS", "PG")):
            a = score_df.loc[score_df["group"] == ref, pop].values
            b = score_df.loc[score_df["group"] == case, pop].values
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append([pop, ref, case, round(t, 3), p, a.mean(), b.mean(),
                         used[pop]])
    res = pd.DataFrame(rows, columns=["population", "ref", "case", "t", "p",
                                      "mean_ref", "mean_case", "n_markers"])
    res["FDR"] = np.nan
    for contrast in res.groupby(["ref", "case"]).groups:
        pass
    for (ref, case), m in res.groupby(["ref", "case"]).groups.items():
        idx = list(m)
        p = res.loc[idx, "p"].values
        n = len(p)
        order = np.argsort(p)
        q = p[order] * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q_inv = np.empty(n)
        q_inv[order] = q
        res.loc[idx, "FDR"] = q_inv
    res.to_csv(os.path.join(TAB, "Table_S7_mcpcounter_stats.csv"),
               index=False, encoding="utf-8")
    print(res.to_string(index=False))

    # Figure S6
    pops = list(score_df.columns[1:])
    order = ["PG", "SS", "NN"]
    means = score_df.groupby("group")[pops].mean().reindex(order)
    fig = plt.figure(figsize=(7.09, 3.7), dpi=600)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.75], wspace=0.5,
                          left=0.13, right=0.97, top=0.90, bottom=0.24)
    axh = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    im = axh.imshow(means.values.T, aspect="auto", cmap="RdBu_r",
                    vmin=-1.2, vmax=1.2)
    axh.set_xticks(range(3))
    axh.set_xticklabels(order, fontsize=6.5)
    axh.set_yticks(range(len(pops)))
    axh.set_yticklabels(pops, fontsize=5.8)
    axh.set_xlabel("Group mean score", fontsize=7)
    cb = fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=5.5)
    axh.set_title("Cell abundance (z)", fontsize=7, pad=2)

    key = ["Neutrophils", "Monocytic lineage", "Myeloid dendritic cells",
           "T cells"]
    palette = {"PG": "#B64342", "SS": "#9A4D8E", "NN": "#767676"}
    xs = np.arange(len(key))
    w = 0.26
    for j, g in enumerate(order):
        vals = [score_df.loc[score_df["group"] == g, p].values for p in key]
        bp = axb.boxplot(vals, positions=xs + (j - 1) * w, widths=w * 0.9,
                         patch_artist=True, showfliers=False,
                         medianprops={"color": "#272727", "lw": 0.8})
        for patch in bp["boxes"]:
            patch.set_facecolor(palette[g])
            patch.set_edgecolor("#272727")
            patch.set_linewidth(0.6)
    axb.set_xticks(xs)
    axb.set_xticklabels([p.replace(" ", "\n") for p in key], fontsize=6)
    axb.set_ylabel("Marker score (log2-normalised)", fontsize=6.5)
    axb.set_title("Key myeloid and lymphoid populations", fontsize=7, pad=2)
    for j, p in enumerate(key):
        for ref, case, yy in (("NN", "PG", 3.3), ("SS", "PG", 3.8)):
            a = score_df.loc[score_df["group"] == ref, p].values
            b = score_df.loc[score_df["group"] == case, p].values
            t, pv = stats.ttest_ind(a, b, equal_var=False)
            if pv < 0.05:
                axb.plot([j - w, j + w], [yy, yy], color="#272727", lw=0.6)
                axb.text(j, yy + 0.12, "*", ha="center", fontsize=6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="s", ls="", color=palette[g],
                      label=g, markersize=6) for g in order]
    axb.legend(handles=handles, loc="upper center", fontsize=6, frameon=False,
               ncol=3, bbox_to_anchor=(0.5, -0.17))
    axh.text(-0.22, 1.06, "a", transform=axh.transAxes, fontsize=9,
             fontweight="bold")
    axb.text(-0.08, 1.06, "b", transform=axb.transAxes, fontsize=9,
             fontweight="bold")
    base = os.path.join(OUT, "Figure_S6_mcpcounter")
    fig.savefig(base + ".svg", bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(base + ".png", dpi=600, bbox_inches="tight")
    fig.canvas.draw()
    fw, fh = fig.canvas.get_width_height()
    issues = []
    for ax in (axh, axb):
        for t in ax.texts:
            bb = t.get_window_extent()
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fw or bb.y1 > fh:
                issues.append(f"OUTSIDE: {t.get_text()[:30]!r}")
    print("QA issues:", len(issues))
    plt.close(fig)


if __name__ == "__main__":
    main()
