# -*- coding: utf-8 -*-
"""Assemble new reproducibility/supplementary tables for v5.1.

Produces:
  S12a GSE280220 key-gene validation
  S12b GSE280220 programme validation
  S13  ODE summary (calibration + uncertainty + graded + combinations)
  S14  Enrichr full raw tables (up/down)
  S15  MR-PRESSO R5 output (CD/UC -> PG)
and mirrors them plus the new analysis scripts into the archived repository.
"""
from __future__ import annotations

import os
import shutil

import pandas as pd

ROOT = os.environ.get("PG_ROOT", r"F:\坏疽性脓皮病\outputs")
TAB = os.path.join(ROOT, "tables")
FIG = os.path.join(ROOT, "figures")
MR = os.path.join(ROOT, "mr")
GSE = os.path.join(MR, "gse280220")
REPO = os.path.join(ROOT, "analysis_repository")
os.makedirs(TAB, exist_ok=True)


def write_csv(df: pd.DataFrame, name: str) -> str:
    path = os.path.join(TAB, name)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def mirror(src: str, subdir: str) -> None:
    if not os.path.exists(src):
        return
    dst = os.path.join(REPO, subdir, os.path.basename(src))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    # ---------- S12a: key genes ----------
    kg = pd.read_csv(os.path.join(GSE, "gse280220_key_genes.csv"))
    keep = ["gene", "PG_mean", "HD_mean", "log2FC_PG_vs_HD", "t", "p"]
    kg["p"] = pd.to_numeric(kg["p"], errors="coerce")
    kg = kg.sort_values("log2FC_PG_vs_HD", ascending=False)
    s12a = write_csv(kg[keep], "Table_S12a_gse280220_key_genes.csv")

    # ---------- S12b: programme validation ----------
    prog = pd.read_csv(os.path.join(GSE, "gse280220_programme_validation.csv"))
    s12b = write_csv(prog, "Table_S12b_gse280220_programmes.csv")

    # ---------- S13: ODE summary ----------
    cal = pd.read_csv(os.path.join(FIG, "pg_network_calibration_ranking.csv"))
    unc = pd.read_csv(os.path.join(FIG, "pg_network_uncertainty_summary.csv"))
    grad = pd.read_csv(os.path.join(FIG, "pg_network_graded_ablation.csv"))
    combo = pd.read_csv(os.path.join(FIG, "pg_network_combo_ablation.csv"))
    unc_piv = unc.pivot(index="intervention", columns="endpoint",
                        values=["median", "p2.5", "p97.5"])
    unc_piv.columns = [f"unc_{a}_{b}" for a, b in unc_piv.columns]
    unc_piv = unc_piv.reset_index()
    ode = cal.merge(unc_piv, left_on="intervention", right_on="intervention",
                    how="left")
    # key graded partial-suppression values (q = 0.3 and 0.7) for NEUT/D
    g_rows = []
    for target in grad["target"].unique():
        row = {"intervention": target}
        for q in (0.3, 0.7):
            sub = grad[(grad["target"] == target) &
                       (grad["remaining_fraction"] == q)]
            if len(sub):
                row[f"q{int(round(q*10))}_NEUT"] = sub.iloc[0]["NEUT"]
                row[f"q{int(round(q*10))}_D"] = sub.iloc[0]["D"]
        g_rows.append(row)
    g30 = pd.DataFrame(g_rows)
    ode = ode.merge(g30, left_on="intervention", right_on="intervention",
                    how="left")
    s13 = write_csv(ode, "Table_S13_ode_ablation_summary.csv")

    # ---------- S14: Enrichr full raw ----------
    tx = os.environ.get("TX_DATA_OUT", r"F:\tx_data\out")
    for label, src in [("up", os.path.join(tx, "enrichr_PG_NN_up.csv")),
                       ("down", os.path.join(tx, "enrichr_PG_NN_down.csv"))]:
        dst = os.path.join(TAB, f"Table_S14_enrichr_full_{label}.csv")
        shutil.copy2(src, dst)

    # ---------- S15: MR-PRESSO R5 ----------
    press = pd.read_csv(os.path.join(MR, "CD_UC_R5_mrpresso.csv"))
    press = press.rename(columns={
        "n_snp": "n_SNP",
        "raw_b": "raw_b",
        "raw_se": "raw_se",
        "raw_t": "raw_t",
        "raw_p": "raw_P",
        "global_RSSobs": "global_RSSobs",
        "global_P": "global_P",
    })
    press.insert(1, "outcome", "PG (FinnGen R5)")
    press["permutations"] = 1000
    press["random_seed"] = 20260814
    s15 = write_csv(press, "Table_S15_mrpresso_R5.csv")

    print("tables:", s12a, s12b, s13, s15, sep="\n  ")
    print("S14 copied:", os.path.exists(os.path.join(TAB, "Table_S14_enrichr_full_up.csv")),
          os.path.exists(os.path.join(TAB, "Table_S14_enrichr_full_down.csv")))

    # ---------- mirror into repository ----------
    for f in ["Table_S12a_gse280220_key_genes.csv",
              "Table_S12b_gse280220_programmes.csv",
              "Table_S13_ode_ablation_summary.csv",
              "Table_S14_enrichr_full_up.csv",
              "Table_S14_enrichr_full_down.csv",
              "Table_S15_mrpresso_R5.csv"]:
        mirror(os.path.join(TAB, f), "results")
    for src in [os.path.join(GSE, "gse280220_gene_validation.csv"),
                os.path.join(GSE, "gse280220_programme_validation.csv"),
                os.path.join(GSE, "gse280220_validation_summary.txt"),
                os.path.join(MR, "CD_UC_R5_mrpresso.txt"),
                os.path.join(MR, "CD_UC_R5_mrpresso.rds"),
                os.path.join(MR, "CD_UC_R5_mrpresso.csv"),
                os.path.join(FIG, "pg_network_calibration_ranking.csv"),
                os.path.join(FIG, "pg_network_calibration_states.csv"),
                os.path.join(FIG, "pg_network_uncertainty_summary.csv"),
                os.path.join(FIG, "pg_network_uncertainty_draws.csv"),
                os.path.join(FIG, "pg_network_graded_ablation.csv"),
                os.path.join(FIG, "pg_network_combo_ablation.csv")]:
        mirror(src, "results")
    for src in [os.path.join(ROOT, "code", "validate_gse280220.py"),
                os.path.join(ROOT, "code", "save_mrpresso_r5.R"),
                os.path.join(ROOT, "code", "make_reproducibility_tables.py")]:
        mirror(src, "code")
    print("repository mirrored")


if __name__ == "__main__":
    main()
