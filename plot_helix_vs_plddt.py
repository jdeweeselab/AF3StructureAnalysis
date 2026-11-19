#!/usr/bin/env python3
"""
plot_helix_vs_plddt.py — summary figure from all_conditions_SS_summary.csv

Top panel:    % of CTD secondary-structure calls that are helix (all models pooled)
Bottom panel: mean pLDDT (± SD) of CTD residues called helix vs. called loop,
              with AlphaFold confidence bands (50, 70, 90)

Usage:  python plot_helix_vs_plddt.py results\\all_conditions_SS_summary.csv [output_basename]
Writes <output_basename>.svg / .pdf / .png  (default: CTD_helix_vs_pLDDT next to the CSV)
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial" if any(
    "Arial" in f.name for f in matplotlib.font_manager.fontManager.ttflist) else "DejaVu Sans"


def nice(r):
    ptm = {"none": "no PTM", "PTM": "PTM"}.get(r["ptm"], r["ptm"])
    dna = {"none": "no DNA", "50+75bp": "50+75 bp"}.get(r["dna"], r["dna"].replace("bp", " bp"))
    return f"{r['protein']} {ptm}, {dna}"


def main():
    src = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_name("CTD_helix_vs_pLDDT")
    R = list(csv.DictReader(open(src, newline="", encoding="utf-8-sig")))
    f = lambda k: np.array([float(r[k]) if r[k] != "" else np.nan for r in R])
    x = np.arange(len(R))
    labels = [nice(r) for r in R]
    # protein group boundaries
    bounds = [i - 0.5 for i in range(1, len(R)) if R[i]["protein"] != R[i - 1]["protein"]]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 5.0), sharex=True)
    a1.bar(x, f("CTD_pct_helix_calls"), color="#d62728", width=0.7)
    a1.set_ylabel("CTD residues called helix\n(% of all calls, 15 models)", fontsize=7)
    a1.set_ylim(0, 100)

    w = 0.38
    ek = dict(lw=0.6, capsize=1.5)
    a2.bar(x - w / 2, f("CTD_pLDDT_helix_mean"), w, yerr=f("CTD_pLDDT_helix_sd"),
           color="#d62728", label="Helix calls", error_kw=ek)
    a2.bar(x + w / 2, f("CTD_pLDDT_loop_mean"), w, yerr=f("CTD_pLDDT_loop_sd"),
           color="#9a9a9a", label="Loop calls", error_kw=ek)
    for y, c in ((90, "#0053D6"), (70, "#65CBF3"), (50, "#FFDB13")):
        a2.axhline(y, color=c, lw=0.8, zorder=0)
    a2.set_ylabel("Mean pLDDT (± SD)", fontsize=7)
    a2.set_ylim(0, 100)
    a2.legend(fontsize=6, frameon=False, ncol=2, loc="upper right")
    a2.set_xticks(x)
    a2.set_xticklabels(labels, rotation=55, ha="right", fontsize=6)

    for a in (a1, a2):
        a.tick_params(labelsize=6, length=2, width=0.5)
        for b in bounds:
            a.axvline(b, color="k", lw=0.6)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_linewidth(0.6)
    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=600, bbox_inches="tight")
    print(f"wrote {out}.svg/.pdf/.png")


if __name__ == "__main__":
    main()
