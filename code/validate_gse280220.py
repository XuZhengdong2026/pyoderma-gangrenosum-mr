# -*- coding: utf-8 -*-
"""Independent validation of the PG lesional hierarchy in GSE280220.

GSE280220 (Di Domizio et al., Nature Communications 2024) profiled the
inflammatory transcriptome of skin lesions with the NanoString nCounter Human
Inflammation panel v2 (GPL19963): 5 pyoderma gangrenosum (PG) lesions,
2 Sweet syndrome lesions and 8 healthy-skin donors.

The script compares PG versus healthy skin for (i) key genes prioritised in
GSE298908, (ii) the seven pre-specified inflammatory programmes, and
(iii) per-gene effect-size correlation with GSE298908.
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.environ.get("PG_ROOT", r"F:\坏疽性脓皮病\outputs")
GEO_DIR = os.path.join(ROOT, "geo")
MR_DIR = os.path.join(ROOT, "mr", "gse280220")
FIG_DIR = os.path.join(ROOT, "figures")
os.makedirs(MR_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

MATRIX = os.path.join(GEO_DIR, "GSE280220_processed.txt.gz")
SKIN_FC = os.path.join(ROOT, "mr", "gut_skin", "skin_log2fc.csv")
PROGRAMME_GENES = os.path.join(ROOT, "mr", "gut_skin", "programme_genes.csv")
PROGRAMME_SUMMARY = os.path.join(ROOT, "mr", "gut_skin", "programme_summary.csv")


def main() -> None:
    mat = pd.read_csv(MATRIX, sep="\t", compression="gzip")
    mat = mat.dropna(subset=["ID_REF"])
    mat = mat[~mat["ID_REF"].str.upper().str.startswith("NEG")]
    val = mat.drop(columns=["GeneID", "Alias"]).groupby("ID_REF", as_index=True).mean()
    val = val.apply(pd.to_numeric, errors="coerce")
    pg_cols = [c for c in val.columns if c.startswith("PG_")]
    hd_cols = [c for c in val.columns if c.startswith("HD_")]
    sw_cols = [c for c in val.columns if c.startswith("SW_")]

    pg = val[pg_cols].astype(float)
    hd = val[hd_cols].astype(float)
    sw = val[sw_cols].astype(float)

    rows = []
    for gene in val.index:
        a, b = pg.loc[gene].values, hd.loc[gene].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        fc = float(np.mean(a) - np.mean(b))
        fc_ss = float(np.mean(a) - np.mean(sw.loc[gene].values))
        rows.append({
            "gene": gene,
            "PG_mean": round(float(np.mean(a)), 4),
            "HD_mean": round(float(np.mean(b)), 4),
            "log2FC_PG_vs_HD": round(fc, 4),
            "t": round(float(t), 4),
            "p": p,
            "log2FC_PG_vs_Sweet": round(fc_ss, 4),
        })
    gene_tab = pd.DataFrame(rows).sort_values("log2FC_PG_vs_HD", ascending=False)
    gene_tab.to_csv(os.path.join(MR_DIR, "gse280220_gene_validation.csv"),
                    index=False, encoding="utf-8-sig")

    key_genes = ["CXCL8", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6",
                 "IL1B", "IL1A", "IL1RN", "CSF3", "CSF3R", "S100A8", "S100A9",
                 "TNF", "IL6", "SOCS3", "JAK2", "STAT3", "NFKB1", "NFKB2",
                 "MMP9", "COL11A1", "IFI44", "IFI6", "IFI44L", "ISG15", "MX1"]
    key = gene_tab[gene_tab["gene"].isin(key_genes)].copy()
    key["available"] = True
    missing = [g for g in key_genes if g not in set(key["gene"])]
    key.to_csv(os.path.join(MR_DIR, "gse280220_key_genes.csv"),
               index=False, encoding="utf-8-sig")

    # ---- programme-level validation ----
    prog_genes = pd.read_csv(PROGRAMME_GENES)
    prog_sum = pd.read_csv(PROGRAMME_SUMMARY)
    prog_rows = []
    for prog, sub in prog_genes.groupby("programme"):
        present = [g for g in sub["gene"] if g in gene_tab.set_index("gene").index]
        if not present:
            continue
        vals = gene_tab.set_index("gene").loc[present, "log2FC_PG_vs_HD"]
        ref = prog_sum.loc[prog_sum["programme"] == prog, "skin_PG_mean"].iloc[0]
        prog_rows.append({
            "programme": prog,
            "n_panel_genes": len(present),
            "n_up_PG": int((vals > 0).sum()),
            "mean_log2FC_GSE280220": round(float(vals.mean()), 4),
            "mean_log2FC_GSE298908": round(float(ref), 4),
            "direction_consistent": bool((vals.mean() > 0) == (ref > 0)),
            "genes": ";".join(present),
        })
    prog_tab = pd.DataFrame(prog_rows)
    prog_tab.to_csv(os.path.join(MR_DIR, "gse280220_programme_validation.csv"),
                    index=False, encoding="utf-8-sig")

    # ---- cross-dataset per-gene correlation ----
    skin = pd.read_csv(SKIN_FC)
    m = gene_tab.merge(skin, left_on="gene", right_on="gene", how="inner")
    m = m[np.isfinite(m["log2FC_PG_vs_HD"]) & np.isfinite(m["fc_PG_vs_NN"])]
    rho, p_rho = stats.spearmanr(m["log2FC_PG_vs_HD"], m["fc_PG_vs_NN"])
    pear, p_pear = stats.pearsonr(m["log2FC_PG_vs_HD"], m["fc_PG_vs_NN"])
    same_sign = float((np.sign(m["log2FC_PG_vs_HD"]) == np.sign(m["fc_PG_vs_NN"])).mean())

    summary = {
        "n_PG": len(pg_cols),
        "n_healthy": len(hd_cols),
        "n_Sweet": len(sw_cols),
        "n_genes_panel": int(len(val)),
        "n_shared_correlation": int(len(m)),
        "spearman_rho": round(float(rho), 4),
        "spearman_p": p_rho,
        "pearson_r": round(float(pear), 4),
        "pearson_p": p_pear,
        "same_sign_proportion": round(same_sign, 4),
        "missing_key_genes": missing,
    }
    with open(os.path.join(MR_DIR, "gse280220_validation_summary.txt"), "w",
              encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}: {v}\n")

    print("=== Key genes (GSE280220 PG vs healthy skin) ===")
    print(key[["gene", "PG_mean", "HD_mean", "log2FC_PG_vs_HD", "p"]].to_string(index=False))
    print("missing:", missing)
    print("\n=== Programme validation ===")
    print(prog_tab.to_string(index=False))
    print("\n=== Cross-dataset correlation ===")
    print(f"n={len(m)} rho={rho:.4f} (p={p_rho:.2e}) pearson={pear:.4f} same_sign={same_sign:.3f}")

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    # panel a: programme means normalised to dataset maximum
    order = prog_tab["programme"].tolist()
    g28 = prog_tab.set_index("programme")["mean_log2FC_GSE298908"].reindex(order)
    g220 = prog_tab.set_index("programme")["mean_log2FC_GSE280220"].reindex(order)
    y = np.arange(len(order))[::-1]
    ax = axes[0]
    ax.barh(y - 0.19, g28.values / g28.abs().max(), height=0.34,
            color="#90A4AE", label="GSE298908 (RNA-seq)")
    ax.barh(y + 0.19, g220.values / g220.abs().max(), height=0.34,
            color="#EF6C00", label="GSE280220 (nCounter)")
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=9)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("Programme mean log2FC (scaled to dataset maximum)")
    ax.set_title("Pre-specified programme direction", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_xlim(-1.02, 1.02)

    ax = axes[1]
    ax.scatter(m["log2FC_PG_vs_HD"], m["fc_PG_vs_NN"], s=12, alpha=0.55,
               color="#1565C0", edgecolors="none")
    ax.axhline(0, color="grey", lw=0.7)
    ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("GSE280220 log2FC (PG vs healthy)")
    ax.set_ylabel("GSE298908 log2FC (PG vs normal)")
    ax.set_title(f"Shared panel genes (n = {len(m)})\nSpearman rho = {rho:.3f}, P = {p_rho:.1e}",
                 fontsize=10)
    fig.suptitle("Independent PG skin transcriptome validation (GSE280220: 5 PG vs 8 healthy skin)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    png = os.path.join(FIG_DIR, "Figure_S8_GSE280220_validation.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "Figure_S8_GSE280220_validation.pdf"),
                bbox_inches="tight")
    from PIL import Image
    im = Image.open(png)
    im.save(os.path.join(FIG_DIR, "Figure_S8_GSE280220_validation.tiff"),
            dpi=(300, 300))
    plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    main()
