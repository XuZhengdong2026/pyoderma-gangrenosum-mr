# -*- coding: utf-8 -*-
"""Figure 5: in silico network model (quantitative grid: time courses + knockout bars)."""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

OUT = r"F:\坏疽性脓皮病\outputs\figures"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, r"F:\坏疽性脓皮病\outputs\code")
import pg_network_knockout as model  # noqa: E402

interventions = {
    "Baseline (PG)": {},
    "CXCL8 KO": {"CXCL8": True},
    "IL1B KO": {"IL1B": True},
    "JAK2 inhibition": {"JAK": 0.1},
    "TNF KO": {"TNF": True},
}
results = {}
for name, ko in interventions.items():
    results[name] = model.simulate(ko)

colors = {
    "Baseline (PG)": "#37474F",
    "CXCL8 KO": "#2E7D32",
    "IL1B KO": "#B71C1C",
    "JAK2 inhibition": "#1565C0",
    "TNF KO": "#6A1B9A",
}

fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), width_ratios=[1.15, 1.0])
ax = axes[0]
t = results["Baseline (PG)"]["sol"].t
for name, c in colors.items():
    sol = results[name]["sol"]
    ax.plot(t, sol.y[model.STATE_NAMES.index("NEUT")], color=c, lw=1.2, label=f"{name}, NEUT")
    ax.plot(t, sol.y[model.STATE_NAMES.index("D")], color=c, lw=1.2, ls="--")
ax.set_xlabel("Time (a.u.)")
ax.set_ylabel("Relative level")
ax.set_title("a  Neutrophils (solid) / ulcer damage (dashed)", fontsize=8, loc="left")
ax.legend(fontsize=5.5, loc="upper left", ncol=1)

ax = axes[1]
pct = pd.read_csv(os.path.join(OUT, "pg_network_knockout_effect.csv"))
order = ["STAT3 KO", "JAK2 inhibition", "CXCL8 KO", "IL1B KO", "TNF KO"]
pct = pct.set_index("intervention").loc[order]
x = np.arange(len(order))
width = 0.36
for j, ep in enumerate(["NEUT", "D"]):
    vals = pct[ep].values
    ax.bar(x + (j - 0.5) * width, vals, width, label=ep,
           color="#2E7D32" if ep == "NEUT" else "#1565C0",
           edgecolor="white", linewidth=0.4)
    for xi, v in zip(x + (j - 0.5) * width, vals):
        ax.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=5)
ax.set_xticks(x)
ax.set_xticklabels(order, fontsize=6)
ax.set_ylabel("Reduction vs baseline (%)")
ax.set_ylim(0, 112)
ax.axhline(0, color="grey", lw=0.7)
ax.set_title("b  In silico knockout effects", fontsize=8, loc="left")
ax.legend(fontsize=6)

fig.tight_layout()
base = os.path.join(OUT, "Figure5_network_model")
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{base}.svg", bbox_inches="tight")
fig.savefig(f"{base}.pdf", bbox_inches="tight")
fig.savefig(f"{base}.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("Figure5 saved", fig.get_size_inches())
