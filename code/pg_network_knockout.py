# -*- coding: utf-8 -*-
"""
PG 炎症网络 in-silico 基因敲除模型
----------------------------------
坏疽性脓皮病（PG）核心炎症回路的确定性 ODE 模型：
组织损伤刺激(S) -> IL-1beta / TNF / NF-kB -> IL-6 -> pSTAT3(JAK2)
-> CXCL8 等中性粒细胞趋化因子 -> 中性粒细胞浸润 -> 组织损伤(溃疡)
SOCS3 为负反馈。

干预（模拟基因敲除/药物）：
  - STAT3 KO          : 消除 pSTAT3 生成
  - JAK2 抑制(模拟baricitinib) : 大幅削弱 IL-6/IL-1beta 对 STAT3 的激活
  - CXCL8 KO          : 消除中性粒细胞趋化信号
  - IL1B KO / TNF KO  : 单细胞因子敲除对照

输出：稳态水平表 + 时间曲线图 + 敲除效果条形图。
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ---------- 参数（相对单位，标定到“PG 样稳态”） ----------
P = {
    "S": 1.0,      # 持续组织损伤/触发刺激
    "k1": 0.50,    # S -> IL1B
    "k2": 0.60,    # NFkB -> IL1B
    "k3": 0.40,    # NFkB -> IL6
    "k4": 0.30,    # pSTAT3 -> IL6
    "k5": 0.40,    # NFkB -> TNF
    "k6": 0.25,    # NEUT -> TNF
    "k7": 0.50,    # IL1B -> NFkB
    "k8": 0.35,    # TNF -> NFkB
    "k9": 0.30,    # NEUT -> NFkB
    "k10": 0.45,   # IL6 -> pSTAT3 (JAK2 依赖)
    "k11": 0.20,   # IL1B -> pSTAT3 (JAK2 依赖)
    "k12": 0.45,   # NFkB -> CXCL8
    "k13": 0.35,   # pSTAT3 -> CXCL8
    "k14": 0.25,   # IL1B -> CXCL8
    "k15": 0.50,   # CXCL8 -> NEUT
    "k16": 0.20,   # TNF -> NEUT
    "k17": 0.35,   # pSTAT3 -> SOCS3
    "k18": 0.15,   # SOCS3 抑制 pSTAT3
    "k19": 0.55,   # NEUT -> 组织损伤 D
    "d1": 0.35, "d2": 0.30, "d3": 0.30, "d4": 0.25,
    "d5": 0.30, "d6": 0.35, "d7": 0.20, "d8": 0.25, "d9": 0.15,
}

STATE_NAMES = ["IL1B", "IL6", "TNF", "NFkB", "pSTAT3", "CXCL8", "NEUT", "SOCS3", "D"]
N0 = len(STATE_NAMES)


def f(x):
    """饱和激活项：x/(1+x)，保证系统有界。"""
    return x / (1.0 + x)


def rhs(t, y, p: dict, ko: dict):
    IL1B, IL6, TNF, NFkB, pSTAT3, CXCL8, NEUT, SOCS3, D = y
    S = p["S"]

    dIL1B = 0.0 if ko.get("IL1B") else p["k1"] * S + p["k2"] * f(NFkB) - p["d1"] * IL1B
    dIL6 = 0.0 if ko.get("IL6") else p["k3"] * f(NFkB) + p["k4"] * f(pSTAT3) - p["d2"] * IL6
    dTNF = 0.0 if ko.get("TNF") else p["k5"] * f(NFkB) + p["k6"] * f(NEUT) - p["d3"] * TNF

    dNFkB = p["k7"] * f(IL1B) + p["k8"] * f(TNF) + p["k9"] * f(NEUT) - p["d4"] * NFkB

    if ko.get("STAT3"):
        dpSTAT3 = -p["d5"] * pSTAT3
    else:
        jak = ko.get("JAK", 1.0)
        dpSTAT3 = jak * (p["k10"] * f(IL6) + p["k11"] * f(IL1B)) - p["k18"] * SOCS3 * pSTAT3 - p["d5"] * pSTAT3

    dCXCL8 = 0.0 if ko.get("CXCL8") else (
        p["k12"] * f(NFkB) + p["k13"] * f(pSTAT3) + p["k14"] * f(IL1B) - p["d6"] * CXCL8
    )
    dNEUT = 0.0 if ko.get("NEUT") else p["k15"] * f(CXCL8) + p["k16"] * f(TNF) - p["d7"] * NEUT

    dSOCS3 = p["k17"] * f(pSTAT3) - p["d8"] * SOCS3
    dD = p["k19"] * f(NEUT) - p["d9"] * D
    return [dIL1B, dIL6, dTNF, dNFkB, dpSTAT3, dCXCL8, dNEUT, dSOCS3, dD]


def simulate(ko: dict, t_end: float = 120.0) -> dict:
    y0 = np.zeros(N0)
    sol = solve_ivp(rhs, (0, t_end), y0, args=(P, ko), method="LSODA",
                    t_eval=np.linspace(0, t_end, 601), rtol=1e-7, atol=1e-9)
    steady = sol.y[:, -1]
    return {"sol": sol, "steady": steady}


def main() -> None:
    out_dir = r"F:\坏疽性脓皮病\outputs\figures"
    os.makedirs(out_dir, exist_ok=True)

    interventions = {
        "Baseline (PG)": {},
        "STAT3 KO": {"STAT3": True},
        "JAK2 inhibition": {"JAK": 0.1},
        "CXCL8 KO": {"CXCL8": True},
        "IL1B KO": {"IL1B": True},
        "TNF KO": {"TNF": True},
    }

    results = {}
    baseline = None
    for name, ko in interventions.items():
        sim = simulate(ko)
        results[name] = sim
        if baseline is None:
            baseline = sim["steady"]

    # ---- 稳态汇总表 ----
    rows = []
    for name, sim in results.items():
        vals = sim["steady"]
        row = {"intervention": name}
        for i, sname in enumerate(STATE_NAMES):
            row[sname] = round(float(vals[i]), 4)
        rows.append(row)
    df = pd.DataFrame(rows)
    csv_path = os.path.join(out_dir, "pg_network_steady_states.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")

    # 相对基线的下降比例（关键终点）
    endpoints = ["IL1B", "TNF", "pSTAT3", "CXCL8", "NEUT", "D"]
    pct = {}
    for name, sim in results.items():
        if name == "Baseline (PG)":
            continue
        vals = sim["steady"]
        base = baseline
        pct[name] = {
            e: round(float(100.0 * (1.0 - vals[STATE_NAMES.index(e)] / base[STATE_NAMES.index(e)])), 1)
            for e in endpoints
        }
    pct_df = pd.DataFrame(pct).T
    pct_csv = os.path.join(out_dir, "pg_network_knockout_effect.csv")
    pct_df.to_csv(pct_csv, index_label="intervention", encoding="utf-8-sig")

    print("=== 稳态水平 ===")
    print(df.to_string(index=False))
    print("\n=== 相对基线下降比例(%) ===")
    print(pct_df.to_string())

    # ---- 图 1：时间曲线 ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = {"Baseline (PG)": "#37474F", "STAT3 KO": "#D32F2F",
              "JAK2 inhibition": "#1976D2", "CXCL8 KO": "#2E7D32"}
    t = results["Baseline (PG)"]["sol"].t
    for name, c in colors.items():
        sim = results[name]
        axes[0].plot(t, sim["sol"].y[STATE_NAMES.index("IL1B")], color=c, lw=1.6, label=f"{name} IL1B")
        axes[0].plot(t, sim["sol"].y[STATE_NAMES.index("pSTAT3")], color=c, lw=1.6, ls="--")
        axes[1].plot(t, sim["sol"].y[STATE_NAMES.index("NEUT")], color=c, lw=1.6, label=f"{name} NEUT")
        axes[1].plot(t, sim["sol"].y[STATE_NAMES.index("D")], color=c, lw=1.6, ls="--")
    axes[0].set_title("IL-1beta (solid) / pSTAT3 (dashed)")
    axes[0].set_xlabel("Time (a.u.)"); axes[0].set_ylabel("Relative level")
    axes[0].legend(fontsize=8)
    axes[1].set_title("Neutrophils (solid) / Ulcer damage (dashed)")
    axes[1].set_xlabel("Time (a.u.)"); axes[1].set_ylabel("Relative level")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "pg_network_timecourse.png"), dpi=300)
    plt.close(fig)

    # ---- 图 2：敲除效果条形图 ----
    fig2, ax = plt.subplots(figsize=(11.5, 6.2))
    x = np.arange(len(pct_df))
    n_ep = len(endpoints)
    width = 0.8 / (n_ep + 0.4)
    offsets = np.linspace(-(n_ep - 1) / 2, (n_ep - 1) / 2, n_ep) * width
    for j, e in enumerate(endpoints):
        ax.bar(x + offsets[j], [pct_df.loc[i, e] for i in pct_df.index],
               width, label=e, edgecolor="white", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(pct_df.index, rotation=18, ha="right", fontsize=9)
    ax.set_ylabel("Reduction vs baseline (%)", fontsize=10)
    ax.set_title("In silico gene knockout effects on the PG inflammatory network",
                 fontsize=11, pad=10)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_ylim(0, 112)
    ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              fontsize=9, frameon=False)
    ax.margins(x=0.02)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, "pg_network_knockout_bars.png"), dpi=300)
    plt.close(fig2)

    # 汇总 JSON
    summary = {
        "model": "PG inflammatory network ODE (IL1B/IL6/TNF/NFkB/pSTAT3/CXCL8/NEUT/SOCS3/D)",
        "steady_states": df.to_dict(orient="records"),
        "percent_reduction_vs_baseline": pct,
        "figures": [
            os.path.join(out_dir, "pg_network_timecourse.png"),
            os.path.join(out_dir, "pg_network_knockout_bars.png"),
        ],
        "tables": [csv_path, pct_csv],
    }
    with open(os.path.join(out_dir, "pg_network_model.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("\nSaved figures/tables to:", out_dir)


if __name__ == "__main__":
    main()
