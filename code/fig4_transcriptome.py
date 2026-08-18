# -*- coding: utf-8 -*-
"""Figure 4: skin transcriptome (quantitative grid: volcano + key-programme dot plots)."""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
DATA = r"F:\tx_data\out"
os.makedirs(OUT, exist_ok=True)

highlight = {
    "CXCL1", "CXCL5", "CXCL6", "CXCL8", "IL1B", "CSF3",
    "S100A8", "S100A9", "COL11A1", "IL1A", "IL36G", "SOCS3",
}

deseq = pd.read_csv(os.path.join(DATA, "deseq_PG_vs_NN.csv"), index_col=0)
deseq = deseq.dropna(subset=["padj"])
deseq["neglog10"] = -np.log10(deseq["padj"].clip(lower=1e-300))

key = pd.read_csv(os.path.join(DATA, "key_genes_mr_crossref.csv"))
pathway_order = [
    "Neutrophil chemokines", "IL-1 axis", "JAK-STAT", "TNF axis",
    "Neutrophil markers", "ECM / fibrosis", "IFN (Sweet contrast)",
]
key["pathway"] = pd.Categorical(key["pathway"], categories=pathway_order, ordered=True)
pathway_colors = {
    "Neutrophil chemokines": "#B71C1C",
    "IL-1 axis": "#E65100",
    "JAK-STAT": "#6A1B9A",
    "TNF axis": "#AD1457",
    "Neutrophil markers": "#1565C0",
    "ECM / fibrosis": "#2E7D32",
    "IFN (Sweet contrast)": "#00838F",
}

fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0],
                      left=0.07, right=0.99, top=0.94, bottom=0.10, wspace=0.28, hspace=0.65)

# Panel a: volcano
ax = fig.add_subplot(gs[0, 0])
up = (deseq["log2FoldChange"] > 1) & (deseq["padj"] < 0.05)
dn = (deseq["log2FoldChange"] < -1) & (deseq["padj"] < 0.05)
ax.scatter(deseq.loc[~up & ~dn, "log2FoldChange"], deseq.loc[~up & ~dn, "neglog10"],
           s=2, color="#BBBBBB", alpha=0.55, rasterized=True)
ax.scatter(deseq.loc[dn, "log2FoldChange"], deseq.loc[dn, "neglog10"],
           s=2, color="#90A4AE", alpha=0.7, rasterized=True)
ax.scatter(deseq.loc[up, "log2FoldChange"], deseq.loc[up, "neglog10"],
           s=3, color="#EF9A9A", alpha=0.7, rasterized=True)
hl = deseq.loc[deseq.index.isin(highlight)]
ax.scatter(hl["log2FoldChange"], hl["neglog10"], s=14, color="#B71C1C",
           edgecolor="black", linewidth=0.3, zorder=5)
for g in ["CXCL1", "CXCL5", "CXCL8", "IL1B", "CSF3", "S100A8", "COL11A1"]:
    r = hl.loc[hl.index == g]
    if len(r):
        ax.annotate(g, (r["log2FoldChange"].iloc[0], r["neglog10"].iloc[0]),
                    fontsize=5.5, xytext=(3, 3), textcoords="offset points", color="#111111")
ax.axhline(-np.log10(0.05), color="#888888", lw=0.6, ls=":")
ax.axvline(1, color="#888888", lw=0.6, ls=":"); ax.axvline(-1, color="#888888", lw=0.6, ls=":")
ax.set_xlabel("log2 fold change (PG vs normal)")
ax.set_ylabel("-log10(FDR)")
ax.set_title("a  DEGs: PG vs normal skin", fontsize=8, loc="left")
ax.text(0.02, 0.98, f"n = {len(deseq):,} genes\n4,285 FDR < 0.05",
        transform=ax.transAxes, fontsize=5.5, va="top", color="#444444")

# Panel b: key programmes PG vs NN
ax = fig.add_subplot(gs[1, 0])
sub = key[key["contrast"] == "PG_vs_NN"].copy()
sub = sub.sort_values(["pathway", "log2FC"], ascending=[True, False])
genes = sub["gene"].tolist()
y = np.arange(len(genes))
for yi, (_, r) in zip(y, sub.iterrows()):
    ax.scatter(r["log2FC"], yi, s=22, color=pathway_colors[r["pathway"]],
               edgecolor="white", linewidth=0.3)
ax.axvline(0, color="#888888", lw=0.7)
ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=5.5)
ax.set_xlim(-4, 11)
ax.set_xlabel("log2FC (PG vs normal)")
ax.set_title("b  Key programmes", fontsize=8, loc="left")

# Panel c: PG vs SS dot plot
ax = fig.add_subplot(gs[1, 1])
sub2 = key[key["contrast"] == "PG_vs_SS"].copy()
sub2 = sub2.sort_values(["pathway", "log2FC"], ascending=[True, False])
genes2 = sub2["gene"].tolist()
y2 = np.arange(len(genes2))
for yi, (_, r) in zip(y2, sub2.iterrows()):
    ax.scatter(r["log2FC"], yi, s=22, color=pathway_colors[r["pathway"]],
               edgecolor="white", linewidth=0.3)
ax.axvline(0, color="#888888", lw=0.7)
ax.set_yticks(y2); ax.set_yticklabels(genes2, fontsize=5.5)
ax.set_xlim(-4, 11)
ax.set_xlabel("log2FC (PG vs Sweet syndrome)")
ax.set_title("c  PG-specific programmes", fontsize=8, loc="left")

# Legend in top-right cell
ax = fig.add_subplot(gs[0, 1])
ax.axis("off")
handles = [plt.Line2D([0], [0], marker="o", ls="", color=c, markersize=7)
           for c in pathway_colors.values()]
ax.legend(handles, pathway_colors.keys(), loc="center", fontsize=6.5,
          title="Gene programme", title_fontsize=7, frameon=False)
ax.text(0.5, 0.78, "Key skin programmes", ha="center", fontsize=8, fontweight="bold")
ax.text(0.5, 0.10, "COL11A1 (ECM) and IFN genes\ndiscriminate PG from Sweet syndrome",
        ha="center", fontsize=6, color="#444444")

fig.savefig(os.path.join(OUT, "Figure4_transcriptome.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure4_transcriptome.svg"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure4_transcriptome.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure4_transcriptome.tiff"), dpi=600,
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
print("Figure4 saved", fig.get_size_inches())
