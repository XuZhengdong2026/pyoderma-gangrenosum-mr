# -*- coding: utf-8 -*-
"""Power for CD -> PG IVW MR based on the observed instrument panel.
Wald-ratio se_i = observed outcome se_i / |beta_exposure_i|."""
from __future__ import annotations

import csv
import math
import os

from scipy import stats

HARM = r"F:\坏疽性脓皮病\outputs\mr\CD_deLange_to_PG_harmonised.csv"
R12 = r"F:\gwas_data\r12_pg_outcome.csv"
OUT = r"F:\坏疽性脓皮病\outputs\tables\Table_S7_mr_power.csv"

OUTCOMES = {
    "FinnGen R5": (280, 208449),
    "FinnGen R12": (703, 470507),
}
ORS = [1.2, 1.3, 1.5, 2.0]


def main():
    rows = []
    with open(HARM, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r.get("mr_keep", "")).strip().lower() in ("true", "1"):
                rows.append(r)
    print("instruments:", len(rows))
    r12 = {}
    with open(R12, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r12.setdefault(r["SNP"], r)
    out_rows = []
    for label, (ncase, nctl) in OUTCOMES.items():
        inv_var = 0.0
        used = 0
        for r in rows:
            try:
                b_e = abs(float(r["beta.exposure"]))
            except Exception:
                continue
            if b_e <= 0:
                continue
            if label == "FinnGen R5":
                try:
                    se_out = float(r["se.outcome"])
                except Exception:
                    continue
            else:
                rec = r12.get(r["SNP"])
                if rec is None:
                    continue
                try:
                    se_out = float(rec["sebeta"])
                except Exception:
                    continue
            if not (se_out > 0):
                continue
            se_wald = se_out / b_e
            inv_var += 1 / (se_wald ** 2)
            used += 1
        print(label, "instruments used:", used)
        se_ivw = math.sqrt(1 / inv_var)
        for or_ in ORS:
            z = math.log(or_) / se_ivw
            power = (1 - stats.norm.cdf(1.96 - z) +
                     stats.norm.cdf(-1.96 - z))
            out_rows.append([label, f"{or_:.1f}", float(round(power * 100, 1)),
                             float(round(se_ivw, 4))])
        # detectable OR at 80% power
        lo, hi = 1.001, 10.0
        for _ in range(80):
            mid = (lo + hi) / 2
            z = math.log(mid) / se_ivw
            p = 1 - stats.norm.cdf(1.96 - z) + stats.norm.cdf(-1.96 - z)
            if p < 0.8:
                lo = mid
            else:
                hi = mid
        out_rows.append([label, "OR at 80% power",
                         float(round((lo + hi) / 2, 2)),
                         float(round(se_ivw, 4))])
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Outcome", "OR", "Power (%)", "IVW SE (log OR)"])
        w.writerows(out_rows)
    for x in out_rows:
        print(x)


if __name__ == "__main__":
    main()
