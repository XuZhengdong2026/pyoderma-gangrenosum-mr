# -*- coding: utf-8 -*-
"""Graded and combinatorial in-silico ablation of the PG inflammatory-network ODE.

Extends the complete-removal (null) ablations with:
1) graded partial suppression (10-90% reduction of each target's production term),
2) pairwise combined ablations evaluated against Bliss independence.

These are mathematical ablations in a stylised model and are NOT equivalent to
experimental gene knockout; they are used for internal-consistency and
hypothesis-prioritisation checks.
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
IDX = {s: i for i, s in enumerate(pgnet.STATE_NAMES)}
ENDPOINTS = ["NEUT", "D"]
OUT_DIR = r"F:\坏疽性脓皮病\outputs\figures"


def reductions(p: dict, ko: dict, base_p: dict | None = None) -> dict:
    base_p = BASE_P if base_p is None else base_p
    pgnet.P = base_p
    base = pgnet.simulate({})["steady"]
    pgnet.P = p
    sim = pgnet.simulate(ko)["steady"]
    return {
        e: float(100.0 * (1.0 - sim[IDX[e]] / base[IDX[e]]))
        for e in ENDPOINTS
    }


def partial_params(target: str, q: float) -> tuple[dict, dict]:
    """Return (params, ko) for remaining production fraction q of a target."""
    p = dict(BASE_P)
    ko: dict = {}
    if target == "CXCL8":
        p["k12"] *= q
        p["k13"] *= q
        p["k14"] *= q
    elif target == "IL1B":
        p["k1"] *= q
        p["k2"] *= q
    elif target == "TNF":
        p["k5"] *= q
        p["k6"] *= q
    elif target == "STAT3":
        p["k10"] *= q
        p["k11"] *= q
    elif target == "JAK2":
        ko = {"JAK": q}
    else:
        raise ValueError(target)
    return p, ko


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    targets = ["CXCL8", "IL1B", "TNF", "STAT3", "JAK2"]
    grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    rows = []
    for t in targets:
        for q in grid:
            p, ko = partial_params(t, q)
            red = reductions(p, ko)
            rows.append({"target": t, "remaining_fraction": q, **red})
    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "pg_network_graded_ablation.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")

    # Pairwise combinations at full ablation
    singles = {t: reductions(*partial_params(t, 0.0)) for t in targets}
    # For combination rows labelled "JAK2i", the clinical-like 90% inhibition
    # (JAK = 0.1) is the relevant single-intervention comparator.
    singles["JAK2"] = reductions(*partial_params("JAK2", 0.1))
    pairs = [
        ("CXCL8+IL1B", "CXCL8", "IL1B"),
        ("CXCL8+TNF", "CXCL8", "TNF"),
        ("IL1B+TNF", "IL1B", "TNF"),
        ("CXCL8+STAT3", "CXCL8", "STAT3"),
        ("IL1B+JAK2i", "IL1B", "JAK2"),
        ("CXCL8+JAK2i", "CXCL8", "JAK2"),
    ]
    combo_rows = []
    for label, a, b in pairs:
        if a == "JAK2":
            pa, koa = dict(BASE_P), {"JAK": 0.1}
        else:
            pa, koa = partial_params(a, 0.0)
        if b == "JAK2":
            koa = {**koa, "JAK": 0.1}
            pb, kob = dict(BASE_P), {"JAK": 0.1}
        else:
            pb, kob = partial_params(b, 0.0)
        # combine parameter scalings (JAK handled through ko)
        pc = dict(BASE_P)
        for k in ("k1", "k2", "k5", "k6", "k10", "k11", "k12", "k13", "k14"):
            pc[k] = pa[k] * pb[k] / BASE_P[k]
        koc = {**koa, **kob}
        obs = reductions(pc, koc)
        row = {"combination": label, "NEUT_observed": obs["NEUT"], "D_observed": obs["D"]}
        for e in ENDPOINTS:
            ea, eb = singles[a][e], singles[b][e]
            row[f"{e}_bliss_expected"] = ea + eb - ea * eb / 100.0
            row[f"{e}_delta_vs_bliss"] = obs[e] - row[f"{e}_bliss_expected"]
        combo_rows.append(row)
    cdf = pd.DataFrame(combo_rows)
    combo_csv = os.path.join(OUT_DIR, "pg_network_combo_ablation.csv")
    cdf.to_csv(combo_csv, index=False, encoding="utf-8-sig")

    # Figure: two dose-response panels + one combination panel
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    colors = {"CXCL8": "#2E7D32", "IL1B": "#C62828", "TNF": "#6A1B9A",
              "STAT3": "#D84315", "JAK2": "#1565C0"}
    for t in targets:
        sub = df[df.target == t].sort_values("remaining_fraction")
        axes[0].plot(sub.remaining_fraction, sub.NEUT, "-o", ms=4, lw=1.8,
                     color=colors[t], label=t)
        axes[1].plot(sub.remaining_fraction, sub.D, "-o", ms=4, lw=1.8,
                     color=colors[t], label=t)
    for ax in (axes[0], axes[1]):
        ax.set_xlabel("Remaining production fraction (q)")
        ax.set_ylabel("Reduction vs baseline (%)")
        ax.legend(fontsize=8, frameon=False)
        ax.axvline(0, color="grey", lw=0.7, ls=":")
        ax.set_xlim(-0.05, 1.05)
    axes[0].set_title("Neutrophil infiltration (NEUT)")
    axes[1].set_title("Ulcer damage (D)")

    x = np.arange(len(cdf))
    w = 0.26
    axes[2].bar(x - w, cdf.NEUT_observed, w, label="Observed NEUT", color="#2E7D32")
    axes[2].bar(x, cdf.NEUT_bliss_expected, w, label="Bliss expected NEUT", color="#A5D6A7")
    axes[2].bar(x + w, cdf.D_observed, w, label="Observed D", color="#1565C0")
    axes[2].bar(x + 2 * w, cdf.D_bliss_expected, w, label="Bliss expected D", color="#90CAF9")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(cdf.combination, rotation=18, ha="right", fontsize=8)
    axes[2].set_ylabel("Reduction vs baseline (%)")
    axes[2].set_title("Combined ablations vs Bliss independence")
    axes[2].legend(fontsize=7.5, frameon=False, loc="upper center",
                   bbox_to_anchor=(0.5, -0.22), ncol=2)
    fig.suptitle("Graded and combinatorial in-silico ablation of the PG inflammatory network",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    fig_path = os.path.join(OUT_DIR, "pg_network_graded_combo.png")
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "graded": df.to_dict(orient="records"),
        "combinations": cdf.to_dict(orient="records"),
        "figure": fig_path,
    }
    with open(os.path.join(OUT_DIR, "pg_network_graded_combo.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("=== Graded ablation: NEUT reduction (%) at q ===")
    print(df.pivot(index="remaining_fraction", columns="target", values="NEUT").round(1).to_string())
    print("\n=== Graded ablation: D reduction (%) at q ===")
    print(df.pivot(index="remaining_fraction", columns="target", values="D").round(1).to_string())
    print("\n=== Combinations (observed vs Bliss) ===")
    print(cdf.to_string(index=False))


if __name__ == "__main__":
    main()
