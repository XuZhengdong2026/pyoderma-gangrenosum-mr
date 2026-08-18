# -*- coding: utf-8 -*-
"""Prepare PG skin (GSE298908) and IBD gut (GSE75214) expression tables."""
import csv
import gzip
import os
import re
from collections import defaultdict

import numpy as np

GEO = r"F:\坏疽性脓皮病\outputs\geo"
OUT = r"F:\坏疽性脓皮病\outputs\mr\gut_skin"
os.makedirs(OUT, exist_ok=True)

# ---------------- PG skin counts (GSE298908) ----------------
def parse_counts(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header  # all columns are sample names (no gene label in header)
        rows = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(header) + 1:
                continue
            gene = p[0]
            vals = np.array([float(x) for x in p[1:]], dtype=float)
            rows.append((gene, vals))
    return samples, rows

print("Parsing GSE298908 counts...")
samples, rows = parse_counts(os.path.join(GEO, "GSE298908_counts_ddsHTSseq_NN_PG_SS.txt.gz"))
def group_of(name):
    if "Pyoderma_Gangrenosum" in name:
        return "PG"
    if "Sweet" in name:
        return "SS"
    if "_NN." in name or ".NN." in name or name.endswith("_NN") or name.endswith(".NN"):
        return "NN"
    return "other"
groups = [group_of(s) for s in samples]
print("n samples:", len(samples), {g: groups.count(g) for g in set(groups)})

mat = np.vstack([r[1] for r in rows])
genes = [r[0] for r in rows]
libsize = mat.sum(axis=0)
cpm = mat / libsize * 1e6
log2cpm = np.log2(cpm + 1)

def mean_log2fc(group_a, group_b):
    a = log2cpm[:, [g == group_a for g in groups]]
    b = log2cpm[:, [g == group_b for g in groups]]
    if a.shape[1] == 0 or b.shape[1] == 0:
        return None
    fc = a.mean(axis=1) - b.mean(axis=1)
    from scipy import stats
    pvals = []
    for i in range(a.shape[0]):
        try:
            t, p = stats.ttest_ind(a[i], b[i], equal_var=False)
            pvals.append(p)
        except Exception:
            pvals.append(np.nan)
    return fc, np.array(pvals)

res_pg = mean_log2fc("PG", "NN")
res_ss = mean_log2fc("SS", "NN")

with open(os.path.join(OUT, "skin_log2fc.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene", "fc_PG_vs_NN", "p_PG_vs_NN", "fc_SS_vs_NN", "p_SS_vs_NN"])
    for i, g in enumerate(genes):
        w.writerow([g,
                    round(float(res_pg[0][i]), 4) if res_pg else "",
                    round(float(res_pg[1][i]), 6) if res_pg else "",
                    round(float(res_ss[0][i]), 4) if res_ss else "",
                    round(float(res_ss[1][i]), 6) if res_ss else ""])
print("skin_log2fc.csv written:", len(genes), "genes")

# ---------------- IBD gut (GSE75214) ----------------
print("Parsing GSE75214 series matrix...")
path = os.path.join(GEO, "GSE75214_series_matrix.txt.gz")
meta = {}
titles = None
table = []
in_table = False
with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            in_table = False
            continue
        if in_table:
            p = line.split("\t")
            table.append(p)
        elif line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!"):
            key = line.split("\t")[0]
            if key not in meta:
                meta[key] = line.split("\t")[1:]

probe_ids = [r[0] for r in table]
if probe_ids and probe_ids[0].strip('"').upper().startswith("ID"):
    table = table[1:]
    probe_ids = [r[0] for r in table]
expr = np.array([[float(x) if x not in ("", "NA") else np.nan for x in r[1:]] for r in table])
print("probes:", len(probe_ids), "samples:", len(titles), "expr shape:", expr.shape)

def cat_of(title):
    t = title.lower()
    if t.startswith("cd_ileum") or t.startswith("cd.ileum"):
        return "CD_ileum"
    if t.startswith("cd_colon"):
        return "CD_colon"
    if t.startswith("uc_colon"):
        return "UC_colon"
    if t.startswith("control_ileum") or t.startswith("controle_ileum"):
        return "CTRL_ileum"
    if t.startswith("control_colon"):
        return "CTRL_colon"
    if t.startswith("control"):
        return "CTRL_other"
    return "other"

cats = [cat_of(t) for t in titles]
print("categories:", {c: cats.count(c) for c in set(cats)})
with open(os.path.join(OUT, "gse75214_samples.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "category"])
    w.writerows(zip(titles, cats))

def gut_fc(cat_a, cat_b):
    ia = [i for i, c in enumerate(cats) if c == cat_a]
    ib = [i for i, c in enumerate(cats) if c == cat_b]
    if not ia or not ib:
        return None
    a = expr[:, ia]
    b = expr[:, ib]
    fc = np.nanmean(a, axis=1) - np.nanmean(b, axis=1)
    return fc

fc_cd_ileum = gut_fc("CD_ileum", "CTRL_ileum")
fc_uc_colon = gut_fc("UC_colon", "CTRL_colon")
fc_cd_colon = gut_fc("CD_colon", "CTRL_colon")

with open(os.path.join(OUT, "gut_probe_fc.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["probe", "fc_CD_ileum", "fc_UC_colon", "fc_CD_colon"])
    for i, p in enumerate(probe_ids):
        w.writerow([p,
                    round(float(fc_cd_ileum[i]), 4) if fc_cd_ileum is not None else "",
                    round(float(fc_uc_colon[i]), 4) if fc_uc_colon is not None else "",
                    round(float(fc_cd_colon[i]), 4) if fc_cd_colon is not None else ""])
print("gut_probe_fc.csv written:", len(probe_ids), "probes")
with open(os.path.join(OUT, "probe_ids.txt"), "w") as f:
    f.write("\n".join(probe_ids))
print("done")
