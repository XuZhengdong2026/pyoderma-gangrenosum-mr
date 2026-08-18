# -*- coding: utf-8 -*-
"""Calibrate the PG inflammatory-network ODE to GSE298908 expression ratios.

The hand-set parameters are re-estimated (log-scale) against the PG-vs-normal
fold changes of the model's measured states, keeping the network structure
fixed. The complete-ablation ranking is then re-evaluated under the
calibrated parameters to test robustness of the CXCL8-first hierarchy.
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SPEC = importlib.util.spec_from_file_location(
    "pgnet", r"F:\坏疽性脓皮病\outputs\code\pg_network_knockout.py"
)
pgnet = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pgnet)

BASE_P = dict(pgnet.P)
STATE_NAMES = pgnet.STATE_NAMES
IDX = {s: i for i, s in enumerate(STATE_NAMES)}
OUT_DIR = r"F:\坏疽性脓皮病\outputs\figures"

# PG vs normal-skin fold changes from GSE298908 (deseq_PG_vs_NN.csv).
# pSTAT3 has no mRNA proxy; STAT3 mRNA is used as a weak anchor.
# NEUT is proxied by the mean log2FC of S100A8, S100A9 and CSF3R.
TARGETS = {
    "IL1B": 2 ** 6.9909,
    "IL6": 2 ** 7.0406,
    "TNF": 2 ** 2.1154,
    "CXCL8": 2 ** 6.7455,
    "SOCS3": 2 ** 4.3692,
    "NEUT": 2 ** ((6.0010 + 6.5084 + 5.9298) / 3.0),
    "NFkB": 2 ** 1.0629,
    "pSTAT3": 2 ** 0.7840,
}
WEIGHTS = {"IL1B": 1.0, "IL6": 1.0, "TNF": 1.0, "CXCL8": 1.0,
           "SOCS3": 1.0, "NEUT": 1.0, "NFkB": 0.5, "pSTAT3": 0.25}
LAMBDA = 0.05


def steady(params_log: np.ndarray) -> np.ndarray:
    p = {k: BASE_P[k] * 10.0 ** float(s) for k, s in zip(BASE_P, params_log)}
    pgnet.P = p
    return pgnet.simulate({})["steady"]


def objective(s):
    st = steady(s)
    errs = []
    for name, target in TARGETS.items():
        val = max(st[IDX[name]], 1e-12)
        errs.append(np.sqrt(WEIGHTS[name]) * (np.log10(val) - np.log10(target)))
    # mild regularisation toward the hand-set parameters
    errs.extend(np.sqrt(LAMBDA) * np.asarray(s))
    return np.asarray(errs)


def reductions(p: dict, ko: dict) -> dict:
    pgnet.P = p
    base = pgnet.simulate({})["steady"]
    sim = pgnet.simulate(ko)["steady"]
    out = {}
    for e in ("NEUT", "D"):
        out[e] = float(100.0 * (1.0 - sim[IDX[e]] / base[IDX[e]]))
    return out


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    x0 = np.zeros(len(BASE_P))
    bounds = (np.full(len(BASE_P), -3.0), np.full(len(BASE_P), 3.0))
    fit = least_squares(objective, x0, bounds=bounds, max_nfev=2000)
    fitted = {k: BASE_P[k] * 10.0 ** float(s) for k, s in zip(BASE_P, fit.x)}

    st0 = steady(np.zeros(len(BASE_P)))
    st1 = steady(fit.x)
    rows = []
    for name in STATE_NAMES:
        rows.append({
            "state": name,
            "target_fold": TARGETS.get(name, np.nan),
            "handset_steady": float(st0[IDX[name]]),
            "calibrated_steady": float(st1[IDX[name]]),
            "calibrated_over_handset": float(st1[IDX[name]] / max(st0[IDX[name]], 1e-12)),
        })
    sdf = pd.DataFrame(rows)
    sdf.to_csv(os.path.join(OUT_DIR, "pg_network_calibration_states.csv"),
               index=False, encoding="utf-8-sig")

    inter = {
        "CXCL8 KO": {"CXCL8": True},
        "IL1B KO": {"IL1B": True},
        "TNF KO": {"TNF": True},
        "STAT3 KO": {"STAT3": True},
        "JAK2 inhibition": {"JAK": 0.1},
    }
    base_rank = {}
    cal_rank = {}
    for name, ko in inter.items():
        base_rank[name] = reductions(BASE_P, ko)
        cal_rank[name] = reductions(fitted, ko)
    rdf = pd.DataFrame({
        "intervention": list(inter),
        "handset_NEUT_pct": [base_rank[k]["NEUT"] for k in inter],
        "handset_D_pct": [base_rank[k]["D"] for k in inter],
        "calibrated_NEUT_pct": [cal_rank[k]["NEUT"] for k in inter],
        "calibrated_D_pct": [cal_rank[k]["D"] for k in inter],
    })
    rdf.to_csv(os.path.join(OUT_DIR, "pg_network_calibration_ranking.csv"),
               index=False, encoding="utf-8-sig")

    summary = {
        "targets": TARGETS,
        "weights": WEIGHTS,
        "fitted_log10_multipliers": {k: float(s) for k, s in zip(BASE_P, fit.x)},
        "states": sdf.to_dict(orient="records"),
        "ranking_handset": base_rank,
        "ranking_calibrated": cal_rank,
    }
    with open(os.path.join(OUT_DIR, "pg_network_calibration.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    print("=== Calibration fit (log10 target vs fitted steady) ===")
    print(sdf.round(4).to_string(index=False))
    print("\n=== Ablation ranking: hand-set vs calibrated ===")
    print(rdf.round(2).to_string(index=False))
    print("\nFitted multipliers (log10):",
          {k: round(float(v), 3) for k, v in zip(BASE_P, fit.x)})


if __name__ == "__main__":
    main()
