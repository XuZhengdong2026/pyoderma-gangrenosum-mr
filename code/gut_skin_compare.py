# -*- coding: utf-8 -*-
"""Compare PG skin (GSE298908) and IBD gut (GSE75214) transcriptomic programmes."""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"F:\坏疽性脓皮病\outputs\mr\gut_skin"

# ---------------- load ----------------
def load_csv(name):
    with open(os.path.join(OUT, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

skin = load_csv("skin_log2fc.csv")
gut = load_csv("gut_probe_fc.csv")
map_rows = load_csv("probe_to_gene.csv")

skin_fc = {r["gene"]: float(r["fc_PG_vs_NN"]) for r in skin if r["fc_PG_vs_NN"] not in ("", "NA")}
probe_gene = {r["PROBEID"]: r["SYMBOL"] for r in map_rows}

gut_by_gene = defaultdict(list)
for r in gut:
    g = probe_gene.get(r["probe"])
    if not g:
        continue
    for col in ("fc_CD_ileum", "fc_UC_colon", "fc_CD_colon"):
        v = r[col]
        if v not in ("", "NA"):
            gut_by_gene[g].append((col, float(v)))

gut_fc = {"CD_ileum": {}, "UC_colon": {}, "CD_colon": {}}
for g, vals in gut_by_gene.items():
    for col, v in vals:
        gut_fc[col.replace("fc_", "")][g] = v

genes = sorted(set(skin_fc) & set(gut_fc["CD_ileum"]) & set(gut_fc["UC_colon"]))
print("shared genes:", len(genes))

sk = np.array([skin_fc[g] for g in genes])
cd = np.array([gut_fc["CD_ileum"][g] for g in genes])
uc = np.array([gut_fc["UC_colon"][g] for g in genes])

rho_cd, p_cd = stats.spearmanr(sk, cd)
rho_uc, p_uc = stats.spearmanr(sk, uc)
print(f"Skin PG vs gut CD ileum: rho={rho_cd:.3f} p={p_cd:.2e}")
print(f"Skin PG vs gut UC colon: rho={rho_uc:.3f} p={p_uc:.2e}")

sign_cd = np.mean(np.sign(sk) == np.sign(cd))
sign_uc = np.mean(np.sign(sk) == np.sign(uc))
print(f"same-sign: CD {sign_cd:.3f}, UC {sign_uc:.3f}")

# ---------------- programmes ----------------
programmes = {
    "Neutrophil chemokines": ["CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8"],
    "IL-1 axis": ["IL1A", "IL1B", "IL1RN", "IL18", "IL36B", "IL36G"],
    "JAK-STAT": ["JAK2", "STAT3", "SOCS3", "IL6", "IL6R", "IL6ST"],
    "TNF/NF-kB": ["TNF", "NFKB1", "NFKB2", "RELA", "TNFAIP3", "TNFAIP6"],
    "Neutrophil/myeloid": ["S100A8", "S100A9", "CSF3R", "FCGR3B", "CD177", "CAMP", "LTF", "MMP9"],
    "Extracellular matrix": ["COL11A1", "COL1A1", "COL3A1", "MMP1"],
    "Interferon-stimulated": ["MX1", "ISG15", "OAS1", "OAS2", "IFIT1", "IFIT3"],
}

prog_rows = []
prog_gene_rows = []
for prog, gs in programmes.items():
    row = {"programme": prog, "n_genes": len(gs)}
    for label, d in (("skin_PG", skin_fc), ("gut_CD_ileum", gut_fc["CD_ileum"]), ("gut_UC_colon", gut_fc["UC_colon"])):
        vals = [d[g] for g in gs if g in d]
        row[label + "_mean"] = float(np.mean(vals)) if vals else np.nan
        row[label + "_n"] = len(vals)
        row[label + "_up"] = sum(1 for v in vals if v > 0)
    prog_rows.append(row)
    for g in gs:
        if g in skin_fc:
            prog_gene_rows.append({
                "programme": prog, "gene": g,
                "skin_PG_fc": skin_fc.get(g, ""),
                "gut_CD_ileum_fc": gut_fc["CD_ileum"].get(g, ""),
                "gut_UC_colon_fc": gut_fc["UC_colon"].get(g, ""),
            })

with open(os.path.join(OUT, "programme_summary.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(prog_rows[0].keys()))
    w.writeheader()
    w.writerows(prog_rows)
with open(os.path.join(OUT, "programme_genes.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(prog_gene_rows[0].keys()))
    w.writeheader()
    w.writerows(prog_gene_rows)

with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as f:
    f.write(f"shared_genes={len(genes)}\n")
    f.write(f"rho_skinPG_vs_gutCD_ileum={rho_cd:.4f}\tp={p_cd:.3e}\n")
    f.write(f"rho_skinPG_vs_gutUC_colon={rho_uc:.4f}\tp={p_uc:.3e}\n")
    f.write(f"same_sign_CD={sign_cd:.4f}\tsame_sign_UC={sign_uc:.4f}\n")

for r in prog_rows:
    print(r["programme"], "skin", round(r["skin_PG_mean"], 2), "CDgut", round(r["gut_CD_ileum_mean"], 2),
          "UCgut", round(r["gut_UC_colon_mean"], 2))

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axes[0]
ax.scatter(cd, sk, s=3, alpha=0.35, color="#1f77b4", edgecolors="none")
ax.set_xlabel("Gut CD ileum log2FC (inflamed vs control)")
ax.set_ylabel("PG skin log2FC (lesional vs normal)")
ax.set_title(f"Shared transcriptomic programmes\nSpearman rho = {rho_cd:.2f}, P = {p_cd:.1e}")
for g in ["CXCL8", "IL1B", "TNF", "S100A8", "S100A9", "CSF3", "COL11A1", "JAK2"]:
    if g in skin_fc and g in gut_fc["CD_ileum"]:
        ax.annotate(g, (gut_fc["CD_ileum"][g], skin_fc[g]), fontsize=7, color="firebrick")
ax.axhline(0, color="grey", lw=0.7); ax.axvline(0, color="grey", lw=0.7)

ax = axes[1]
labels = list(prog_rows)
names = [r["programme"] for r in prog_rows]
x = np.arange(len(names))
width = 0.28
skin_mean = [r["skin_PG_mean"] for r in prog_rows]
cd_mean = [r["gut_CD_ileum_mean"] for r in prog_rows]
uc_mean = [r["gut_UC_colon_mean"] for r in prog_rows]
ax.bar(x - width, skin_mean, width, label="PG skin", color="#d62728")
ax.bar(x, cd_mean, width, label="Gut CD ileum", color="#1f77b4")
ax.bar(x + width, uc_mean, width, label="Gut UC colon", color="#ff7f0e")
ax.set_xticks(x)
ax.set_xticklabels([n.replace(" ", "\n") for n in names], fontsize=7)
ax.set_ylabel("Mean log2FC")
ax.axhline(0, color="grey", lw=0.7)
ax.legend(fontsize=7)
ax.set_title("Pre-specified gene programmes")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure_gut_skin_shared.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Figure_gut_skin_shared.tiff"), dpi=300)
print("figure written")
