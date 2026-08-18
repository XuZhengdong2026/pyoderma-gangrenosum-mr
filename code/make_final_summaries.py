# -*- coding: utf-8 -*-
"""Regenerate final reproducible summary tables (MR summary, S01 combined, S15 combined)."""
import csv
import math
import os
import shutil

OUT = r"F:\坏疽性脓皮病\outputs"
MR = os.path.join(OUT, "mr")
TAB = os.path.join(OUT, "tables")
REPO = os.path.join(OUT, "analysis_repository", "results")


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def or_ci(b, se):
    return math.exp(b), math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)


# ---------- 1. Final MR results summary (strict LD-clumped, R5 + R12) ----------
main_sources = [
    ("IBD", "R5", os.path.join(MR, "IBD_ldclump_mr_results.csv")),
    ("CD", "R5", os.path.join(MR, "CD_deLange_to_PG_mr_results.csv")),
    ("UC", "R5", os.path.join(MR, "UC_deLange_to_PG_mr_results.csv")),
    ("IBD", "R12", os.path.join(MR, "IBD_R12_mr_results.csv")),
    ("CD", "R12", os.path.join(MR, "CD_deLange_R12_mr_results.csv")),
    ("UC", "R12", os.path.join(MR, "UC_deLange_R12_mr_results.csv")),
]
summary = []
for exposure, freeze, path in main_sources:
    for r in read_csv(path):
        method = r["method"].strip()
        b = float(r["b"])
        se = float(r["se"])
        p = float(r["pval"])
        orv, lo, hi = or_ci(b, se)
        summary.append({
            "exposure": exposure,
            "data_freeze": freeze,
            "method": method,
            "nsnp": r["nsnp"],
            "OR": round(orv, 4),
            "OR_lo": round(lo, 4),
            "OR_hi": round(hi, 4),
            "beta": round(b, 6),
            "se": round(se, 6),
            "pval": "%.4g" % p,
        })
write_csv(os.path.join(TAB, "MR_results_summary.csv"), summary,
          ["exposure", "data_freeze", "method", "nsnp", "OR", "OR_lo", "OR_hi",
           "beta", "se", "pval"])

# ---------- 2. Instrument F statistics ----------
f_sources = [
    ("IBD", "R5", os.path.join(MR, "IBD_ldclump_harmonised.csv")),
    ("CD", "R5", os.path.join(MR, "CD_deLange_to_PG_harmonised.csv")),
    ("UC", "R5", os.path.join(MR, "UC_deLange_to_PG_harmonised.csv")),
    ("IBD", "R12", os.path.join(MR, "IBD_R12_harmonised.csv")),
    ("CD", "R12", os.path.join(MR, "CD_deLange_R12_harmonised.csv")),
    ("UC", "R12", os.path.join(MR, "UC_deLange_R12_harmonised.csv")),
]
f_rows = []
for exposure, freeze, path in f_sources:
    fs = []
    for r in read_csv(path):
        if str(r.get("mr_keep", "")).strip().upper() != "TRUE":
            continue
        try:
            b = float(r["beta.exposure"])
            se = float(r["se.exposure"])
        except Exception:
            continue
        if se > 0:
            fs.append((b / se) ** 2)
    f_rows.append({
        "exposure": exposure,
        "data_freeze": freeze,
        "n_instruments": len(fs),
        "mean_F": round(sum(fs) / len(fs), 2) if fs else None,
        "min_F": round(min(fs), 2) if fs else None,
    })
write_csv(os.path.join(TAB, "MR_instrument_F_stats.csv"), f_rows,
          ["exposure", "data_freeze", "n_instruments", "mean_F", "min_F"])
print("F stats:", [(r["exposure"], r["data_freeze"], r["n_instruments"],
                    r["mean_F"], r["min_F"]) for r in f_rows])

# ---------- 3. Supplementary Table S1 (R12 forward + relaxed reverse, Steiger-filtered) ----------
s1 = []
for r in read_csv(os.path.join(MR, "R12_sensitivity_summary.csv")):
    orv = float(r["OR"])
    lo = float(r["lo"])
    hi = float(r["hi"])
    b = math.log(orv)
    se = (math.log(hi) - math.log(lo)) / 3.92
    s1.append({
        "analysis": "Forward MR (FinnGen R12)",
        "exposure": r["exposure"],
        "outcome": "PG",
        "threshold": "",
        "method": "IVW",
        "nsnp": r["nsnp"],
        "beta": round(b, 6),
        "se": round(se, 6),
        "pval": r["ivw_p"],
        "OR": r["OR"],
        "OR_lo": r["lo"],
        "OR_hi": r["hi"],
        "Q": r["Q"],
        "Q_p": r["Q_p"],
        "egger_int_p": r["egger_int_p"],
    })
rev_dir = os.path.join(MR, "reverse_relaxed")
rev_files = [
    ("IBD", "PG_to_IBD_p5e-05_steiger_mr.csv", "5e-05"),
    ("CD", "PG_to_CD_p5e-05_steiger_mr.csv", "5e-05"),
    ("UC", "PG_to_UC_p5e-05_steiger_mr.csv", "5e-05"),
    ("CD", "PG_to_CD_p5e-06_steiger_mr.csv", "5e-06"),
    ("UC", "PG_to_UC_p5e-06_steiger_mr.csv", "5e-06"),
]
for exposure, fname, thr in rev_files:
    path = os.path.join(rev_dir, fname)
    if not os.path.exists(path):
        continue
    for r in read_csv(path):
        b = float(r["b"])
        se = float(r["se"])
        p = float(r["pval"])
        orv, lo, hi = or_ci(b, se)
        s1.append({
            "analysis": "Reverse MR (PG to exposure, Steiger-filtered)",
            "exposure": exposure,
            "outcome": exposure,
            "threshold": thr,
            "method": r["method"].strip(),
            "nsnp": r["nsnp"],
            "beta": round(b, 6),
            "se": round(se, 6),
            "pval": "%.4g" % p,
            "OR": round(orv, 4),
            "OR_lo": round(lo, 4),
            "OR_hi": round(hi, 4),
            "Q": "",
            "Q_p": "",
            "egger_int_p": "",
        })
write_csv(os.path.join(TAB, "S01_R12_sensitivity_reverse_MR.csv"), s1,
          ["analysis", "exposure", "outcome", "threshold", "method", "nsnp",
           "beta", "se", "pval", "OR", "OR_lo", "OR_hi", "Q", "Q_p",
           "egger_int_p"])
print("S1 rows:", len(s1))

# ---------- 4. Supplementary Table S15 (R5 + R12 MR-PRESSO, CD and UC) ----------
s15 = []
for r in read_csv(os.path.join(TAB, "Table_S15_mrpresso_R5.csv")):
    s15.append({
        "exposure": r["exposure"],
        "outcome": "PG (FinnGen R5)",
        "n_SNP": r["n_SNP"],
        "raw_b": r["raw_b"],
        "raw_se": r["raw_se"],
        "raw_t": r["raw_t"],
        "raw_P": r["raw_P"],
        "global_RSSobs": r["global_RSSobs"],
        "global_P": r["global_P"],
        "outliers": r["outliers"],
        "permutations": r["permutations"],
        "random_seed": r["random_seed"],
    })
r12_rows = [
    ("CD", 76, 0.1115871, 0.04996063, 2.2335, 0.02849758, 97.39611, 0.065),
    ("UC", 51, 0.01195645, 0.0517261, 0.2311492, 0.8181422, 44.5252, 0.768),
]
for exposure, n, b, se, t, p, rss, gp in r12_rows:
    s15.append({
        "exposure": exposure,
        "outcome": "PG (FinnGen R12)",
        "n_SNP": n,
        "raw_b": b,
        "raw_se": se,
        "raw_t": t,
        "raw_P": p,
        "global_RSSobs": rss,
        "global_P": gp,
        "outliers": "none detected",
        "permutations": 1000,
        "random_seed": "not recorded",
    })
write_csv(os.path.join(TAB, "Table_S15_mrpresso.csv"), s15,
          ["exposure", "outcome", "n_SNP", "raw_b", "raw_se", "raw_t", "raw_P",
           "global_RSSobs", "global_P", "outliers", "permutations",
           "random_seed"])
print("S15 rows:", len(s15))

# ---------- 5. Copy new reproducibility files into the archived repository ----------
os.makedirs(REPO, exist_ok=True)
for name in ["MR_results_summary.csv", "MR_instrument_F_stats.csv",
             "S01_R12_sensitivity_reverse_MR.csv", "Table_S15_mrpresso.csv"]:
    shutil.copy2(os.path.join(TAB, name), os.path.join(REPO, name))
for fname in os.listdir(rev_dir):
    if fname.endswith("_steiger_mr.csv"):
        shutil.copy2(os.path.join(rev_dir, fname), os.path.join(REPO, fname))
print("REPO_COPY_DONE")
