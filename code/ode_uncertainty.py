# -*- coding: utf-8 -*-
"""Parameter-uncertainty analysis for the PG inflammatory-network ODE.

Perturbs every rate/decay constant by +/-30% (log-uniform) across 800 draws,
reruns baseline and all in-silico interventions, and reports the median and
95% interval of the reduction in neutrophils (NEUT) and ulcer damage (D).
"""
from __future__ import annotations

import importlib.util
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SPEC = importlib.util.spec_from_file_location(
    "pgnet", r"F:\坏疽性脓皮病\outputs\code\pg_network_knockout.py"
)
pgnet = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pgnet)

BASE_P = dict(pgnet.P)
INTERVENTIONS = {
    "CXCL8 KO": {"CXCL8": True},
    "IL1B KO": {"IL1B": True},
    "TNF KO": {"TNF": True},
    "STAT3 KO": {"STAT3": True},
    "JAK2 inhibition": {"JAK": 0.1},
}
ENDPOINTS = ["NEUT", "D"]
IDX = {s: i for i, s in enumerate(pgnet.STATE_NAMES)}


def reductions(p: dict, ko: dict) -> dict:
    pgnet.P = p
    base = pgnet.simulate({})["steady"]
    sim = pgnet.simulate(ko)["steady"]
    out = {}
    for e in ENDPOINTS:
        out[e] = 100.0 * (1.0 - sim[IDX[e]] / base[IDX[e]])
    return out


def deterministic_reference():
    return {
        name: reductions(BASE_P, ko)
        for name, ko in INTERVENTIONS.items()
    }


def main() -> None:
    rng = np.random.default_rng(20260814)
    n = 800
    params = list(BASE_P.keys())
    log_fac = np.log(1.3)
    draws = []
    for rep in range(n):
        factors = np.exp(rng.uniform(-log_fac, log_fac, size=len(params)))
        p = {k: BASE_P[k] * float(f) for k, f in zip(params, factors)}
        row = {"draw": rep + 1}
        ok = True
        for name, ko in INTERVENTIONS.items():
            try:
                r = reductions(p, ko)
                for e in ENDPOINTS:
                    row[f"{name}|{e}"] = r[e]
                # baseline stability check
                pgnet.P = p
                st = pgnet.simulate({})["steady"]
                if not np.all(np.isfinite(st)) or np.any(st < 0):
                    ok = False
            except Exception:
                ok = False
                break
        if ok:
            draws.append(row)
    df = pd.DataFrame(draws)
    print("valid draws:", len(df))

    ref = deterministic_reference()
    summary_rows = []
    for name, ko in INTERVENTIONS.items():
        for e in ENDPOINTS:
            col = f"{name}|{e}"
            s = df[col]
            summary_rows.append({
                "intervention": name,
                "endpoint": e,
                "deterministic": round(ref[name][e], 2),
                "median": round(float(s.median()), 2),
                "p2.5": round(float(s.quantile(0.025)), 2),
                "p97.5": round(float(s.quantile(0.975)), 2),
                "min": round(float(s.min()), 2),
                "max": round(float(s.max()), 2),
            })
    summary = pd.DataFrame(summary_rows)
    print(summary.to_string(index=False))

    out_dir = r"F:\坏疽性脓皮病\outputs\figures"
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, "pg_network_uncertainty_draws.csv"),
              index=False, encoding="utf-8-sig")
    summary.to_csv(os.path.join(out_dir, "pg_network_uncertainty_summary.csv"),
                   index=False, encoding="utf-8-sig")

    # Forest-style plot: median with 95% interval, deterministic as diamond
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    for ax, e in zip(axes, ENDPOINTS):
        sub = summary[summary["endpoint"] == e].set_index("intervention")
        order = list(INTERVENTIONS)
        y = np.arange(len(order))[::-1]
        med = [sub.loc[i, "median"] for i in order]
        lo = [sub.loc[i, "p2.5"] for i in order]
        hi = [sub.loc[i, "p97.5"] for i in order]
        det = [sub.loc[i, "deterministic"] for i in order]
        ax.errorbar(med, y, xerr=[np.array(med) - np.array(lo),
                                  np.array(hi) - np.array(med)],
                    fmt="o", color="#1565C0", ms=6, capsize=4, lw=1.4,
                    label="Median (95% interval)")
        ax.plot(det, y, marker="D", ls="", color="#C62828", ms=7,
                label="Deterministic reference")
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=9)
        ax.set_xlabel(f"{e} reduction vs baseline (%)")
        ax.set_title(e, fontsize=11)
        ax.legend(fontsize=8, loc="lower right", frameon=False)
    fig.suptitle("ODE parameter-uncertainty analysis (±30% on all rate/decay constants, n=800)",
                 fontsize=12, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(out_dir, "pg_network_uncertainty.png"), dpi=300)
    plt.close(fig)

    with open(os.path.join(out_dir, "pg_network_uncertainty.json"), "w",
              encoding="utf-8") as f:
        json.dump({
            "n_draws": len(df),
            "perturbation": "+/-30% log-uniform on k1-k19, d1-d9",
            "summary": summary.to_dict(orient="records"),
        }, f, ensure_ascii=False, indent=2)
    print("\nSaved uncertainty summary/figures to", out_dir)


if __name__ == "__main__":
    main()
