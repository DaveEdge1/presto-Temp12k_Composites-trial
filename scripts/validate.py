#!/usr/bin/env python3
"""Validation for the Temperature 12k composites (holocene_da-style).

Compares the reconstructed per-method + consensus GMST against the published
Temperature 12k curves (Kaufman et al. 2020), bundled under
reference_data/published/<method>_published.csv (age_bp, median, lo, hi).
Reports Pearson r and the coefficient of efficiency (CE) per method (computed on
full-record-centred anomalies so the reference-period convention cancels) and
the mid-Holocene magnitude vs the paper's Table 1, and writes a self-contained
HTML report modeled on holocene_da's validate_holocene_da.py:

  results/validation/index.html          report page (surfaced on GitHub Pages)
  results/validation/gmst_validation.png  overlay of mine vs published
  results/validation/comparison.json      metrics
  results/figures/validation.png          same overlay (static-viz / figures)
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# paper Table 1 mid-Holocene (6.5-5.5 ka) GMST vs 1800-1900 (median, degC)
TABLE1 = {"scc": 0.50, "dcc": 0.50, "gam": 0.44, "cps": 1.08, "paico": 0.42}
METHODS = ["scc", "dcc", "gam", "cps", "paico"]
METHOD_LABEL = {"scc": "SCC", "dcc": "DCC", "gam": "GAM", "cps": "CPS",
                "paico": "PaiCo"}


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


def n_overlap(a, b):
    return int((np.isfinite(a) & np.isfinite(b)).sum())


def recenter(x):
    return x - np.nanmean(x)


def _fmt(v, fmt="{:.3f}", dash="—"):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return dash
    try:
        return fmt.format(float(v))
    except (TypeError, ValueError):
        return str(v)


def build_html(results, age_range, n_members):
    """Render the holocene_da-style validation report."""
    cons = results.get("consensus", {})

    # ── metric cards: consensus first, then per-method r ──
    cards = []
    if "pearson_r" in cons:
        cards.append((_fmt(cons["pearson_r"], "{:.3f}"), "Consensus GMST R vs published"))
    if "coef_efficiency" in cons:
        cards.append((_fmt(cons["coef_efficiency"], "{:.3f}"), "Consensus GMST CE vs published"))
    if "midHolocene_mine" in cons:
        cards.append((_fmt(cons["midHolocene_mine"], "{:.3f}") + " °C",
                      "Consensus mid-Holocene ΔT (6.5–5.5 ka)"))
    for m in METHODS:
        if m in results and "pearson_r" in results[m]:
            cards.append((_fmt(results[m]["pearson_r"], "{:.3f}"),
                          f"{METHOD_LABEL[m]} GMST R vs published"))
    cards_html = "\n".join(
        f'''    <div class="metric-card">
      <div class="value">{val}</div>
      <div class="label">{label}</div>
    </div>''' for val, label in cards)

    # ── per-method + consensus table ──
    rows = []
    order = METHODS + ["consensus"]
    for m in order:
        if m not in results:
            continue
        rec = results[m]
        label = "Consensus" if m == "consensus" else METHOD_LABEL[m]
        published = rec.get("midHolocene_published")
        if m != "consensus":
            published = rec.get("midHolocene_published", TABLE1.get(m))
            pub_note = "" if "midHolocene_published" in rec else " *"
        else:
            pub_note = ""
        ce = rec.get("coef_efficiency")
        ce_color = "#16a34a" if (ce is not None and np.isfinite(ce) and ce > 0.5) \
            else ("#d97706" if (ce is not None and np.isfinite(ce) and ce > 0) else "#dc2626")
        rows.append(
            f'''    <tr>
      <td>{label}</td>
      <td>{rec.get("n_points", "—")}</td>
      <td>{_fmt(rec.get("midHolocene_mine"))} °C</td>
      <td>{_fmt(published)}{pub_note} °C</td>
      <td>{_fmt(rec.get("pearson_r"), "{:.4f}")}</td>
      <td style="color:{ce_color}; font-weight:600;">{_fmt(rec.get("coef_efficiency"), "{:.4f}")}</td>
    </tr>''')
    table_rows = "\n".join(rows) or (
        '    <tr><td colspan="6" style="text-align:center;color:#6b7280;">'
        'No published reference curves found in reference_data/published/.</td></tr>')

    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Temperature 12k Validation</title>
  <style>
    :root {{ --accent: #4682b4; --bg: #f7f8fa; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           max-width: 1100px; margin: 0 auto; padding: 24px; color: #1a1a1a;
           background: var(--bg); }}
    h1 {{ border-bottom: 3px solid var(--accent); padding-bottom: 12px; font-size: 1.8rem; }}
    h2 {{ color: #374151; margin-top: 36px; font-size: 1.3rem;
          border-left: 4px solid var(--accent); padding-left: 12px; }}
    p {{ line-height: 1.6; color: #4b5563; }}
    table {{ border-collapse: collapse; margin: 16px 0; width: 100%;
             background: white; border-radius: 8px; overflow: hidden;
             box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
    th, td {{ border: 1px solid #e5e7eb; padding: 10px 16px; text-align: left; }}
    th {{ background: #f3f4f6; font-weight: 600; font-size: 0.9rem;
          text-transform: uppercase; letter-spacing: 0.03em; color: #6b7280; }}
    img {{ max-width: 100%; margin: 12px 0; border: 1px solid #e5e7eb;
           border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }}
    .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                    gap: 16px; margin: 16px 0; }}
    .metric-card {{ background: white; padding: 20px; border-radius: 8px;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.08); text-align: center; }}
    .metric-card .value {{ font-size: 2rem; font-weight: 700; color: var(--accent); }}
    .metric-card .label {{ font-size: 0.85rem; color: #6b7280; margin-top: 4px; }}
    .back {{ margin-top: 32px; }}
    code {{ background: #eef1f4; padding: 1px 6px; border-radius: 3px; font-size: 0.9em; }}
    .note {{ font-size: 0.85rem; color: #6b7280; }}
  </style>
</head>
<body>
  <h1>Temperature 12k Validation Report</h1>
  <p>Validation of this multi-method composite reconstruction against the
     <strong>published Kaufman et al. (2020) Temperature 12k</strong> curves
     (<a href="https://doi.org/10.1038/s41597-020-0530-7">doi:10.1038/s41597-020-0530-7</a>),
     bundled from the NOAA archive under <code>reference_data/published/</code>.
     Pearson correlation (R) and the Nash–Sutcliffe coefficient of efficiency (CE)
     are computed on full-record-centred anomalies (so the reference-period
     convention cancels) of the ensemble-median GMST. CE = 1 is perfect;
     CE = 0 equals climatology; CE &lt; 0 is worse than climatology.</p>

  <div class="metric-grid">
{cards_html}
  </div>

  <h2>GMST Validation Metrics</h2>
  <p>Per-method and consensus global-mean surface temperature against the
     published Temp12k curves over their common age range
     ({age_range[0]:,}–{age_range[1]:,} yr BP). The mid-Holocene column is the
     mean anomaly over 6.5–5.5 ka.</p>
  <table>
    <tr><th>Method</th><th>n pts</th><th>mid-Holocene (mine)</th>
        <th>mid-Holocene (published)</th><th>Pearson R</th><th>CE</th></tr>
{table_rows}
  </table>
  <p class="note">* where no published mid-Holocene value is bundled, the
     paper's Table 1 value is shown for reference.</p>

  <h2>GMST Time Series: Reconstruction vs Published</h2>
  <p>Per-method (coloured) and consensus (black) reconstruction medians (solid)
     overlaid on the published Temp12k curves (dotted). The shaded band is the
     consensus 5–95% ensemble spread. X-axis runs from the oldest age (left)
     to present (right).</p>
  <img src="gmst_validation.png" alt="GMST reconstruction vs published">

  <p class="note">Pooled multi-method ensemble: {n_members:,} members.</p>
  <p class="back"><a href="../index.html">&larr; Back to results</a></p>
</body>
</html>"""


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
    vdir = out / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    results = {}
    fig, ax = plt.subplots(figsize=(11, 6))
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
        ax.plot(kyr, mine, color=colors[m], lw=1.4, label=f"{METHOD_LABEL[m]} (mine)")
        ref_p = refdir / f"{m}_published.csv"
        if ref_p.exists():
            ref = pd.read_csv(ref_p)
            ref_med = np.interp(age, ref["age_bp"].to_numpy(float),
                                ref["median"].to_numpy(float))
            rec["pearson_r"] = round(pearson_r(recenter(mine), recenter(ref_med)), 4)
            rec["coef_efficiency"] = round(coef_efficiency(recenter(ref_med), recenter(mine)), 4)
            rec["n_points"] = n_overlap(mine, ref_med)
            ax.plot(kyr, ref_med, color=colors[m], lw=1.0, ls=":", alpha=0.8,
                    label=f"{METHOD_LABEL[m]} (published)")
        results[m] = rec

    n_members = int(d["n_members"].iloc[0]) if "n_members" in d.columns else 0
    if "mean" in d.columns:
        cons = d["mean"].to_numpy(float)
        if {"lo_95", "hi_95"}.issubset(d.columns):
            ax.fill_between(kyr, d["lo_95"].to_numpy(float), d["hi_95"].to_numpy(float),
                            color="0.5", alpha=0.18, label="consensus 5–95%")
        ax.plot(kyr, cons, color="k", lw=2.4, label="consensus (mine)")
        rec = {"midHolocene_mine": round(float(np.nanmean(cons[mh])), 3)}
        cons_p = refdir / "consensus_published.csv"
        if cons_p.exists():
            cp = pd.read_csv(cons_p)
            cref = np.interp(age, cp["age_bp"].to_numpy(float), cp["median"].to_numpy(float))
            rec["midHolocene_published"] = round(float(np.nanmean(cref[mh])), 3)
            rec["pearson_r"] = round(pearson_r(recenter(cons), recenter(cref)), 4)
            rec["coef_efficiency"] = round(coef_efficiency(recenter(cref), recenter(cons)), 4)
            rec["n_points"] = n_overlap(cons, cref)
            ax.plot(kyr, cref, color="k", lw=1.2, ls=":", alpha=0.8, label="consensus (published)")
        results["consensus"] = rec

    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))
    ax.set_xlabel("Age (ka BP)"); ax.set_ylabel("Temperature anomaly (°C)")
    ax.set_title("Temperature 12k: reconstruction (solid) vs published (dotted)")
    ax.legend(fontsize=7, ncol=3, loc="lower center")
    ax.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(vdir / "gmst_validation.png", dpi=130)
    shutil.copyfile(vdir / "gmst_validation.png", out / "figures" / "validation.png")
    plt.close(fig)

    (vdir / "comparison.json").write_text(json.dumps(results, indent=2))

    age_range = (int(np.nanmin(age)), int(np.nanmax(age)))
    (vdir / "index.html").write_text(build_html(results, age_range, n_members),
                                     encoding="utf-8")

    print("[validate] metrics:")
    for k, v in results.items():
        print(f"  {k}: {v}")
    print(f"[validate] wrote {vdir/'index.html'} + gmst_validation.png + comparison.json")


if __name__ == "__main__":
    main()
