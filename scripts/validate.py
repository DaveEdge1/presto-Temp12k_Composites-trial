#!/usr/bin/env python3
"""Validation for the Temperature 12k composites (holocene_da-style).

Compares the reconstructed per-method + consensus GMST against the published
Temperature 12k curves (Kaufman et al. 2020), bundled under
reference_data/published/<method>_published.csv (age_bp, median, lo, hi).
Reports Pearson r and the coefficient of efficiency (CE) per method (computed on
full-record-centred anomalies so the reference-period convention cancels) and
the mid-Holocene magnitude vs the paper's Table 1, and writes:

  results/validation/comparison.json   metrics
  results/figures/validation.png       overlay of mine vs published

Modeled on holocene_da's validate_holocene_da.py (which likewise validates
against the published Temp12k reconstruction).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# paper Table 1 mid-Holocene (6.5-5.5 ka) GMST vs 1800-1900 (median, degC)
TABLE1 = {"scc": 0.50, "dcc": 0.50, "gam": 0.44, "cps": 1.08, "paico": 0.42}
METHODS = ["scc", "dcc", "gam", "cps", "paico"]


def pearson_r(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() >= 5 else float("nan")


def coef_efficiency(obs, pred):
    m = np.isfinite(obs) & np.isfinite(pred)
    if m.sum() < 5:
        return float("nan")
    o, p = obs[m], pred[m]
    ss_tot = np.sum((o - np.mean(o)) ** 2)
    return float(1.0 - np.sum((o - p) ** 2) / ss_tot) if ss_tot else float("nan")


def recenter(x):
    return x - np.nanmean(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recon", required=True, help="results/reconstruction.csv")
    ap.add_argument("--refdir", required=True, help="reference_data/published")
    ap.add_argument("--out-dir", required=True, help="results dir")
    args = ap.parse_args()

    d = pd.read_csv(args.recon)
    age = d["age_bp"].to_numpy(float)
    mh = (age >= 5500) & (age <= 6500)
    refdir = Path(args.refdir)
    out = Path(args.out_dir)
    (out / "validation").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    results = {}
    fig, ax = plt.subplots(figsize=(10, 5))
    kyr = age / 1000.0
    colors = {"scc": "tab:purple", "dcc": "tab:orange", "gam": "tab:green",
              "cps": "tab:blue", "paico": "tab:red"}

    for m in METHODS:
        col = f"{m}_median"
        if col not in d.columns:
            continue
        mine = d[col].to_numpy(float)
        rec = {"midHolocene_mine": round(float(np.nanmean(mine[mh])), 3),
               "table1": TABLE1[m]}
        ax.plot(kyr, mine, color=colors[m], lw=1.4, label=f"{m.upper()} (mine)")
        ref_p = refdir / f"{m}_published.csv"
        if ref_p.exists():
            ref = pd.read_csv(ref_p)
            ref_med = np.interp(age, ref["age_bp"].to_numpy(float),
                                ref["median"].to_numpy(float))
            rec["pearson_r"] = round(pearson_r(recenter(mine), recenter(ref_med)), 4)
            rec["coef_efficiency"] = round(coef_efficiency(recenter(ref_med), recenter(mine)), 4)
            ax.plot(kyr, ref_med, color=colors[m], lw=1.0, ls=":", alpha=0.8,
                    label=f"{m.upper()} (published)")
        results[m] = rec

    if "mean" in d.columns:
        cons = d["mean"].to_numpy(float)
        ax.plot(kyr, cons, color="k", lw=2.4, label="consensus (mine)")
        rec = {"midHolocene_mine": round(float(np.nanmean(cons[mh])), 3)}
        cons_p = refdir / "consensus_published.csv"
        if cons_p.exists():
            cp = pd.read_csv(cons_p)
            cref = np.interp(age, cp["age_bp"].to_numpy(float), cp["median"].to_numpy(float))
            rec["midHolocene_published"] = round(float(np.nanmean(cref[mh])), 3)
            rec["pearson_r"] = round(pearson_r(recenter(cons), recenter(cref)), 4)
            rec["coef_efficiency"] = round(coef_efficiency(recenter(cref), recenter(cons)), 4)
            ax.plot(kyr, cref, color="k", lw=1.2, ls=":", alpha=0.8, label="consensus (published)")
        results["consensus"] = rec

    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))
    ax.set_xlabel("Age (ka BP)"); ax.set_ylabel("Temperature anomaly (°C)")
    ax.set_title("Temperature 12k: reconstruction (solid) vs published (dotted)")
    ax.legend(fontsize=7, ncol=3, loc="lower center")
    ax.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(out / "figures" / "validation.png", dpi=130)
    plt.close(fig)

    (out / "validation" / "comparison.json").write_text(json.dumps(results, indent=2))
    print("[validate] metrics:")
    for k, v in results.items():
        print(f"  {k}: {v}")
    print(f"[validate] wrote {out/'validation'/'comparison.json'} + figures/validation.png")


if __name__ == "__main__":
    main()
