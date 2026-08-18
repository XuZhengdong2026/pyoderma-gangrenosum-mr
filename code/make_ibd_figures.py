# -*- coding: utf-8 -*-
"""Build Figure 2 (MR forest) and Figure 3 (transcriptome montage) for IBD v3.0."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG = r"F:\坏疽性脓皮病\outputs\figures"
MR = r"F:\坏疽性脓皮病\outputs\mr"

# ---------------- Figure 2 ----------------
main = [
    ("IBD -> PG (R5)", 1.16, 0.98, 1.37),
    ("CD -> PG (R5)", 1.24, 1.05, 1.46),
    ("UC -> PG (R5)", 0.98, 0.81, 1.18),
    ("IBD -> PG (R12)", 1.12, 1.00, 1.24),
    ("CD -> PG (R12)", 1.12, 1.01, 1.23),
    ("UC -> PG (R12)", 1.01, 0.91, 1.13),
]

spec_rows = list(csv.DictReader(open(os.path.join(MR, "skin_specificity_mr.csv"), encoding="utf-8")))
spec = []
for r in spec_rows:
    if r["method"] != "Inverse variance weighted":
        continue
    spec.append((f"{r['exposure']} -> {r['outcome'].replace(' (European GWAS)', '').replace(' (FinnGen R5)', '')}",
                 float(r["or"]), float(r["or_lci"]), float(r["or_uci"])))

fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
for ax, data, title in ((axes[0], main, "Main MR: IBD/CD/UC on PG"),
                        (axes[1], spec, "Specificity: CD/UC on other skin outcomes")):
    data = data[::-1]
    y = np.arange(len(data))
    for i, (lab, orv, lci, uci) in enumerate(data):
        col = "#d62728" if (lci > 1 or uci < 1) else "#1f77b4"
        ax.errorbar(orv, i, xerr=[[orv - lci], [uci - orv]], fmt="s", color=col,
                    markersize=5, linewidth=1.4, capsize=2)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([d[0] for d in data], fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(0.55, 2.4)
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_title(title, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Figure2_IBD_MR_forest.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Figure2_IBD_MR_forest.tiff"), dpi=300)
print("Figure 2 written")

# ---------------- Figure 3 montage ----------------
from PIL import Image
top = Image.open(os.path.join(FIG, "Figure4_transcriptome.png"))
bottom = Image.open(os.path.join(MR, "gut_skin", "Figure_gut_skin_shared.png"))
w = max(top.width, bottom.width)
top2 = top.resize((w, int(top.height * w / top.width)))
bot2 = bottom.resize((w, int(bottom.height * w / bottom.width)))
canvas = Image.new("RGB", (w, top2.height + bot2.height + 20), "white")
canvas.paste(top2, (0, 0))
canvas.paste(bot2, (0, top2.height + 20))
canvas.save(os.path.join(FIG, "Figure3_transcriptome_gutskin.png"))
canvas.save(os.path.join(FIG, "Figure3_transcriptome_gutskin.tiff"))
print("Figure 3 written", canvas.size)
