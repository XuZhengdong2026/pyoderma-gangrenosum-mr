# -*- coding: utf-8 -*-
"""PROGENy pathway-activity analysis on GSE298908 (PG / Sweet / normal skin)
using the official human model (top 100 genes per pathway by p-value)."""
from __future__ import annotations

import csv
import math
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
MODEL = r"F:\gwas_data\progeny_model_human.csv"
OUT = r"F:\坏疽性脓皮病\outputs\figures"
TAB = r"F:\坏疽性脓皮病\outputs\tables"


def read_counts(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        samples = header  # file has no gene-name header; first field is a sample
        genes = []
        rows = []
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


def deseq_size_factors(counts):
    """Median-of-ratios normalisation (DESeq2-style), returns norm + 1 (log2)."""
    gm = np.exp(np.mean(np.log(counts + 0.5), axis=1))
    gm[gm <= 0] = np.nan
    ratio = counts / gm[:, None]
    with np.errstate(invalid="ignore"):
        sf = np.nanmedian(np.where(counts > 0, ratio, np.nan), axis=0)
    sf = np.where(np.isfinite(sf) & (sf > 0), sf, np.nanmedian(sf))
    norm = counts / sf
    # +1 pseudocount keeps the log transform defined for zero counts
    return np.log2(norm + 1)


def build_model(path, top=100):
    df = pd.read_csv(path)
    df = df[["gene", "pathway", "weight", "p.value"]]
    df["gene"] = df["gene"].str.upper()
    genes = sorted(set(df["gene"]))
    pathways = sorted(set(df["pathway"]))
    keep = []
    for pth in pathways:
        sub = df[df["pathway"] == pth].sort_values("p.value").head(top)
        keep.append(sub)
    sub = pd.concat(keep)
    mat = pd.pivot_table(sub, index="gene", columns="pathway",
                         values="weight", aggfunc="first", fill_value=0.0)
    mat = mat.reindex(index=genes, columns=pathways, fill_value=0.0)
    return mat


def main():
    counts, genes, samples = read_counts(COUNTS)
    expr = deseq_size_factors(counts)
    groups = np.array([group_of(s) for s in samples])
    model = build_model(MODEL)
    common = sorted(set(genes) & set(model.index))
    idx_e = [genes.index(g) for g in common]
    M = model.loc[common].values
    scores = expr[idx_e, :].T @ M  # samples x pathways
    scores = (scores - scores.mean(axis=0)) / scores.std(axis=0, ddof=1)
    score_df = pd.DataFrame(scores, index=samples, columns=model.columns)
    score_df.insert(0, "group", groups)
    score_df.to_csv(os.path.join(TAB, "Table_S4_progeny_scores.csv"),
                    encoding="utf-8")

    rows = []
    for pth in model.columns:
        for ref, case in (("NN", "PG"), ("SS", "PG")):
            a = score_df.loc[score_df["group"] == ref, pth].values
            b = score_df.loc[score_df["group"] == case, pth].values
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append([pth, f"{ref} vs {case}", round(t, 3),
                         f"{p:.3g}", a.mean(), b.mean()])
    res = pd.DataFrame(rows, columns=["pathway", "contrast", "t", "p",
                                      "mean_ref", "mean_case"])
    # BH FDR within contrasts
    res["FDR"] = np.nan
    for contrast in res["contrast"].unique():
        m = res["contrast"] == contrast
        p = res.loc[m, "p"].astype(float).values
        n = len(p)
        order = np.argsort(p)
        p_sorted = p[order]
        q = p_sorted * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q_inv = np.empty(n)
        q_inv[order] = q
        res.loc[m, "FDR"] = q_inv
    res.to_csv(os.path.join(TAB, "Table_S5_progeny_stats.csv"),
               index=False, encoding="utf-8")
    print(res.sort_values(["contrast", "FDR"]).to_string(index=False))
    print("common genes:", len(common))

    # ---- Figure S5 ----
    fig = plt.figure(figsize=(7.09, 3.9), dpi=600)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.55], wspace=0.42,
                          left=0.09, right=0.97, top=0.90, bottom=0.16)
    axh = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    order = ["PG", "SS", "NN"]
    means = score_df.groupby("group")[model.columns].mean().reindex(order)
    im = axh.imshow(means.values.T, aspect="auto", cmap="RdBu_r",
                    vmin=-1.5, vmax=1.5)
    axh.set_xticks(range(3))
    axh.set_xticklabels(["PG", "SS", "NN"], fontsize=6.5)
    axh.set_yticks(range(len(model.columns)))
    axh.set_yticklabels(model.columns, fontsize=6)
    axh.set_xlabel("Group mean pathway activity", fontsize=7)
    cb = fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=5.5)
    axh.set_title("Pathway activity (z)", fontsize=7, pad=2)

    key = ["NFkB", "TNFa", "JAK-STAT", "PI3K"]
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
    axb.set_xticklabels(key, fontsize=7)
    axb.set_ylabel("Pathway activity (z)", fontsize=7)
    axb.set_title("Key model-relevant pathways", fontsize=7, pad=2)
    axb.tick_params(axis="x", length=0)
    for j, p in enumerate(key):
        for ref, case, yy in (("NN", "PG", 3.6), ("SS", "PG", 4.1)):
            a = score_df.loc[score_df["group"] == ref, p].values
            b = score_df.loc[score_df["group"] == case, p].values
            t, pv = stats.ttest_ind(a, b, equal_var=False)
            if pv < 0.05:
                axb.plot([j - w, j + w], [yy, yy], color="#272727", lw=0.6)
                axb.text(j, yy + 0.15, "*" if pv < 0.05 else "",
                         ha="center", fontsize=6)
    axb.set_ylim(-2.5, 4.8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="s", ls="", color=palette[g],
                      label=g, markersize=6) for g in order]
    axb.legend(handles=handles, loc="upper right", fontsize=6, frameon=False)
    axh.text(-0.16, 1.06, "a", transform=axh.transAxes, fontsize=9,
             fontweight="bold")
    axb.text(-0.12, 1.06, "b", transform=axb.transAxes, fontsize=9,
             fontweight="bold")
    base = os.path.join(OUT, "Figure_S5_progeny")
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
