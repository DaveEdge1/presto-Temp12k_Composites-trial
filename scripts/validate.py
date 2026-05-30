#!/usr/bin/env python3
"""Validation for the Temperature 12k composites.

Compares the reconstructed per-method + consensus GMST against the published
Kaufman et al. (2020) Temperature 12k curves, bundled under
reference_data/published/<method>_published.csv (age_bp, median, lo, hi).

Both my curves and the published curves are anchored to 0 at 100 yr BP (the
official NOAA archival convention: per ensemble member remove the full 12 ka
mean, then set each method's median to 0 at 100 yr BP / 1800-1900 CE).
Metrics are computed on absolute (not recentered) anomalies so amplitude and
offset errors are visible.

Outputs:
  results/validation/index.html         report page (surfaced on GitHub Pages)
  results/validation/comparison.json    full metric dict per method
  results/validation/gmst_validation.png  combined overlay (mine vs published, all methods)
  results/validation/method_<m>.png     per-method overlay (median + 5-95% band)
  results/figures/validation.png        same combined overlay for static viz
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Paper Table 1 mid-Holocene (5.5-6.5 ka) GMST relative to 1800-1900 (median, degC)
TABLE1 = {"scc": 0.50, "dcc": 0.50, "gam": 0.44, "cps": 1.08, "paico": 0.42,
          "consensus": 0.70}
METHODS = ["scc", "dcc", "gam", "cps", "paico"]
METHOD_LABEL = {"scc": "SCC", "dcc": "DCC", "gam": "GAM", "cps": "CPS",
                "paico": "PaiCo"}
COLORS = {"scc": "tab:purple", "dcc": "tab:orange", "gam": "tab:green",
          "cps": "tab:blue", "paico": "tab:red", "consensus": "black"}


def anchor(age, y, at=100.0):
    o = np.argsort(age)
    return y - np.interp(at, age[o], y[o])


def window_mean(age, y, lo, hi):
    m = (age >= lo) & (age <= hi)
    return float(np.nanmean(y[m])) if m.any() else float("nan")


def flat_run(y):
    best = run = 1
    for i in range(1, len(y)):
        run = run + 1 if y[i] == y[i - 1] else 1
        best = max(best, run)
    return best


def my_curve(method, recon_csv, methods_dir):
    """Return (age, median, q05, q95) for a method.

    For per-method ensembles read methods/<m>_global.csv if present
    (binAges + ensemble columns); else fall back to the consensus reconstruction
    CSV's <m>_median column with no spread.
    """
    g = Path(methods_dir) / f"{method}_global.csv"
    if g.exists():
        df = pd.read_csv(g)
        age = df["binAges"].to_numpy(float)
        ens = df.drop(columns=["binAges"]).to_numpy(float)
        o = np.argsort(age)
        return (age[o], np.nanmedian(ens, 1)[o],
                np.nanpercentile(ens, 5, 1)[o],
                np.nanpercentile(ens, 95, 1)[o])
    d = pd.read_csv(recon_csv)
    age = d["age_bp"].to_numpy(float)
    o = np.argsort(age); age = age[o]
    col = f"{method}_median"
    if col not in d.columns:
        return None
    return age, d[col].to_numpy(float)[o], None, None


def consensus_curve(recon_csv):
    d = pd.read_csv(recon_csv)
    age = d["age_bp"].to_numpy(float); o = np.argsort(age); age = age[o]
    lo = d["lo_95"].to_numpy(float)[o] if "lo_95" in d.columns else None
    hi = d["hi_95"].to_numpy(float)[o] if "hi_95" in d.columns else None
    return age, d["mean"].to_numpy(float)[o], lo, hi


def pub_curve(method, refdir):
    p = Path(refdir) / f"{method}_published.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    age = df["age_bp"].to_numpy(float); o = np.argsort(age)
    return (age[o], df["median"].to_numpy(float)[o],
            df["lo"].to_numpy(float)[o], df["hi"].to_numpy(float)[o])


def compute_metrics(age, med, q05, q95, pub, pub_lo, pub_hi):
    """All metrics anchored at 100 BP, absolute."""
    pmed = np.interp(age, pub[0][np.argsort(pub[0])], pub[1][np.argsort(pub[0])])
    plo = np.interp(age, pub[0][np.argsort(pub[0])], pub_lo[np.argsort(pub[0])])
    phi = np.interp(age, pub[0][np.argsort(pub[0])], pub_hi[np.argsort(pub[0])])
    me = anchor(age, med); pu = anchor(age, pmed)
    msk = np.isfinite(me) & np.isfinite(pu)
    if msk.sum() < 5:
        return None
    r = float(np.corrcoef(me[msk], pu[msk])[0, 1])
    bias = float(np.nanmean(me - pu))
    rmse = float(np.sqrt(np.nanmean((me - pu) ** 2)))
    maxd = float(np.nanmax(np.abs(me - pu)))
    amp = float(np.nanstd(me) / np.nanstd(pu)) if np.nanstd(pu) > 0 else float("nan")
    mh_m = window_mean(age, me, 5500, 6500)
    mh_p = window_mean(age, pu, 5500, 6500)
    c12_m = window_mean(age, me, 11500, 12000)
    c12_p = window_mean(age, pu, 11500, 12000)
    spread = float("nan")
    if q05 is not None and q95 is not None:
        s_m = float(np.nanmean(q95 - q05))
        s_p = float(np.nanmean(phi - plo))
        spread = s_m / s_p if s_p > 0 else float("nan")
    return dict(r=r, bias=bias, RMSE=rmse, maxD=maxd, amp=amp,
                midHol_mine=mh_m, midHol_pub=mh_p, twelve_ka_mine=c12_m, twelve_ka_pub=c12_p,
                spread=spread, flat=int(flat_run(med)), n_points=int(msk.sum()))


def plot_method(method, age, med, q05, q95, pub, pub_lo, pub_hi, out_path):
    o_age = age
    o_pub_age = pub[0]
    me = anchor(o_age, med); pu = anchor(o_pub_age, pub[1])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    kyr = o_age / 1000.0; pkyr = o_pub_age / 1000.0
    col = COLORS.get(method, "black")
    if q05 is not None and q95 is not None:
        ax.fill_between(kyr, anchor(o_age, q05), anchor(o_age, q95),
                        color=col, alpha=0.18, label="mine 5–95%")
    ax.plot(kyr, me, color=col, lw=2.2, label=f"{method.upper() if method != 'consensus' else 'Consensus'} (mine)")
    ax.fill_between(pkyr, anchor(o_pub_age, pub_lo), anchor(o_pub_age, pub_hi),
                    color="0.4", alpha=0.15, label="published 5–95%")
    ax.plot(pkyr, pu, color="0.2", lw=1.4, ls=":", label="published median")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))
    ax.set_xlabel("Age (ka BP)"); ax.set_ylabel("Temperature anomaly (°C)")
    ax.set_title(f"{method.upper() if method != 'consensus' else 'Consensus'}: reconstruction vs published, anchored at 100 BP")
    ax.legend(fontsize=8, loc="lower center", ncol=2)
    ax.grid(alpha=0.3, lw=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_combined(curves_with_pub, out_path):
    """Single combined overlay: all per-method medians + consensus.

    curves_with_pub[name] = (age, med, q05, q95, pub_age, pub_med, pub_lo, pub_hi)
    """
    fig, ax = plt.subplots(figsize=(11, 6))
    for name, (age, med, q05, q95, pub_age, pub_med, pub_lo, pub_hi) in curves_with_pub.items():
        kyr = age / 1000.0; pkyr = pub_age / 1000.0
        col = COLORS.get(name, "black")
        lw = 2.2 if name == "consensus" else 1.3
        ax.plot(kyr, anchor(age, med), color=col, lw=lw,
                label=f"{name.upper() if name != 'consensus' else 'Consensus'} (mine)")
        ax.plot(pkyr, anchor(pub_age, pub_med), color=col, lw=lw - 0.4, ls=":", alpha=0.7,
                label=f"{name.upper() if name != 'consensus' else 'Consensus'} (pub)")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlim(12, 0)
    ax.set_xlabel("Age (ka BP)"); ax.set_ylabel("Temperature anomaly (°C)")
    ax.set_title("Temperature 12k: reconstruction (solid) vs published (dotted), anchored at 100 BP")
    ax.legend(fontsize=7, ncol=3, loc="lower center")
    ax.grid(alpha=0.3, lw=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _fmt(v, fmt="{:.3f}", dash="—"):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return dash
    try:
        return fmt.format(float(v))
    except (TypeError, ValueError):
        return str(v)


def build_html(results, age_range, n_members, per_method_pngs):
    cons = results.get("consensus", {})
    cards = []
    if "r" in cons:        cards.append((_fmt(cons["r"], "{:.3f}"), "Consensus R vs published"))
    if "maxD" in cons:     cards.append((_fmt(cons["maxD"]), "Consensus maxΔ (°C)"))
    if "amp" in cons:      cards.append((_fmt(cons["amp"], "{:.2f}"), "Consensus amplitude (mine/pub)"))
    if "spread" in cons and np.isfinite(cons.get("spread", float("nan"))):
        cards.append((_fmt(cons["spread"], "{:.2f}"), "Consensus envelope ratio"))
    if "midHol_mine" in cons:
        cards.append((_fmt(cons["midHol_mine"]) + " °C", "Consensus mid-Holocene (mine)"))
    if "midHol_pub" in cons:
        cards.append((_fmt(cons["midHol_pub"]) + " °C", "Consensus mid-Holocene (published)"))
    cards_html = "\n".join(
        f'''    <div class="metric-card">
      <div class="value">{val}</div>
      <div class="label">{label}</div>
    </div>''' for val, label in cards)

    rows = []
    for m in METHODS + ["consensus"]:
        if m not in results: continue
        rec = results[m]
        label = "Consensus" if m == "consensus" else METHOD_LABEL[m]
        spread_v = rec.get("spread")
        rows.append(
            f'''    <tr>
      <td><strong>{label}</strong></td>
      <td>{_fmt(rec.get("r"), "{:.3f}")}</td>
      <td>{_fmt(rec.get("bias"), "{:+.3f}")}</td>
      <td>{_fmt(rec.get("RMSE"))}</td>
      <td>{_fmt(rec.get("maxD"))}</td>
      <td>{_fmt(rec.get("amp"), "{:.2f}")}</td>
      <td>{_fmt(rec.get("midHol_mine"), "{:.2f}")} | {_fmt(rec.get("midHol_pub"), "{:.2f}")} | {_fmt(TABLE1.get(m), "{:.2f}")}</td>
      <td>{_fmt(rec.get("twelve_ka_mine"), "{:.2f}")} | {_fmt(rec.get("twelve_ka_pub"), "{:.2f}")}</td>
      <td>{_fmt(spread_v, "{:.2f}") if (spread_v is not None and np.isfinite(spread_v)) else "—"}</td>
      <td>{rec.get("flat", "—")}</td>
    </tr>''')
    table_rows = "\n".join(rows) or '    <tr><td colspan="10" style="text-align:center;color:#6b7280;">No published reference curves.</td></tr>'

    method_imgs = "\n".join(
        f'  <h3>{METHOD_LABEL.get(m, m).upper() if m != "consensus" else "Consensus"}</h3>\n'
        f'  <img src="{png.name}" alt="{m} overlay">'
        for m, png in per_method_pngs.items())

    return f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Temperature 12k Validation</title>
<style>
  :root {{ --accent: #4682b4; --bg: #f7f8fa; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1180px; margin: 0 auto; padding: 24px; color: #1a1a1a; background: var(--bg); }}
  h1 {{ border-bottom: 3px solid var(--accent); padding-bottom: 12px; font-size: 1.8rem; }}
  h2 {{ color: #374151; margin-top: 36px; font-size: 1.3rem; border-left: 4px solid var(--accent); padding-left: 12px; }}
  h3 {{ color: #374151; margin-top: 24px; font-size: 1.05rem; }}
  p  {{ line-height: 1.6; color: #4b5563; }}
  table {{ border-collapse: collapse; margin: 16px 0; width: 100%; background: white;
           border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.08);
           font-size: 0.92rem; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 8px 12px; text-align: left; }}
  th {{ background: #f3f4f6; font-weight: 600; font-size: 0.78rem;
        text-transform: uppercase; letter-spacing: 0.03em; color: #6b7280; }}
  img {{ max-width: 100%; margin: 8px 0 20px; border: 1px solid #e5e7eb;
         border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }}
  .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                  gap: 16px; margin: 16px 0; }}
  .metric-card {{ background: white; padding: 16px; border-radius: 8px;
                  box-shadow: 0 1px 3px rgba(0,0,0,0.08); text-align: center; }}
  .metric-card .value {{ font-size: 1.6rem; font-weight: 700; color: var(--accent); }}
  .metric-card .label {{ font-size: 0.8rem; color: #6b7280; margin-top: 4px; }}
  code {{ background: #eef1f4; padding: 1px 6px; border-radius: 3px; font-size: 0.9em; }}
  .note {{ font-size: 0.85rem; color: #6b7280; }}
  .back {{ margin-top: 32px; }}
</style>
</head>
<body>
<h1>Temperature 12k Validation Report</h1>
<p>Validation of this multi-method composite reconstruction against the
   <strong>published Kaufman et al. (2020) Temperature 12k</strong> curves
   (<a href="https://doi.org/10.1038/s41597-020-0530-7">doi:10.1038/s41597-020-0530-7</a>),
   bundled from the NOAA archive under <code>reference_data/published/</code>.
   Both curves are anchored to 0 at 100 yr BP (NOAA archival convention) and
   metrics computed on <em>absolute</em> anomalies so amplitude and offset
   errors are visible — unlike r/CE on recentered series, which hide them.</p>

<div class="metric-grid">
{cards_html}
</div>

<h2>Validation metrics</h2>
<p>Per-method and consensus comparison against the published curves. Both
   anchored to 0 at 100 yr BP. <code>n pts</code> ages 0–12 ka common to both.
   <strong>r</strong>: Pearson correlation of anchored medians.
   <strong>bias / RMSE / maxΔ</strong>: of (mine − published), °C.
   <strong>amp</strong>: std(mine)/std(pub) — 1.0 = matched amplitude.
   <strong>midHol</strong>: mean over 5.5–6.5 ka — mine | published | paper Table 1.
   <strong>12 ka</strong>: mean over 11.5–12 ka — mine | published.
   <strong>spread</strong>: mean(q95−q05) mine/published — 1.0 = matched uncertainty.
   <strong>flat</strong>: longest run of identical consecutive median values (clipping flag).</p>

<table>
  <tr><th>Method</th><th>r</th><th>bias</th><th>RMSE</th><th>maxΔ</th><th>amp</th>
      <th>midHol (m | p | T1)</th><th>12 ka (m | p)</th><th>spread</th><th>flat</th></tr>
{table_rows}
</table>

<h2>Per-method overlays</h2>
<p>Reconstruction (coloured line + shaded 5–95% spread) vs published curve
   (dotted + grey shaded). X-axis: oldest age left, present right.</p>
{method_imgs}

<h2>Combined overlay</h2>
<img src="gmst_validation.png" alt="Combined GMST overlay">

<p class="note">Pooled multi-method ensemble: {n_members:,} members. Age range
   in the comparison: {age_range[0]:,}–{age_range[1]:,} yr BP.</p>
<p class="back"><a href="../index.html">&larr; Back to results</a></p>
</body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recon", required=True, help="results/reconstruction.csv")
    ap.add_argument("--refdir", required=True, help="reference_data/published")
    ap.add_argument("--out-dir", required=True, help="results dir")
    args = ap.parse_args()

    out = Path(args.out_dir); vdir = out / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    methods_dir = Path(args.recon).parent / "methods"
    refdir = Path(args.refdir)

    # curves[name] = (age, med, q05, q95, pub_age, pub_med, pub_lo, pub_hi)
    curves = {}
    for m in METHODS:
        cu = my_curve(m, args.recon, methods_dir)
        pub = pub_curve(m, refdir)
        if cu is None or pub is None:
            continue
        curves[m] = (*cu, *pub)
    age, cmed, clo, chi = consensus_curve(args.recon)
    pub_c = pub_curve("consensus", refdir)
    if pub_c is not None:
        curves["consensus"] = (age, cmed, clo, chi, *pub_c)

    results = {}
    per_method_pngs = {}
    for name, (a, med, q05, q95, pa, pmed, plo, phi) in curves.items():
        mtr = compute_metrics(a, med, q05, q95, (pa, pmed, plo, phi), plo, phi)
        if mtr is None:
            continue
        results[name] = mtr
        png = vdir / f"method_{name}.png"
        plot_method(name, a, med, q05, q95, (pa, pmed, plo, phi), plo, phi, png)
        per_method_pngs[name] = png

    plot_combined(curves, vdir / "gmst_validation.png")
    shutil.copyfile(vdir / "gmst_validation.png", out / "figures" / "validation.png")

    (vdir / "comparison.json").write_text(json.dumps(results, indent=2))

    d = pd.read_csv(args.recon)
    n_members = int(d["n_members"].iloc[0]) if "n_members" in d.columns else 0
    age_all = d["age_bp"].to_numpy(float)
    age_range = (int(np.nanmin(age_all)), int(np.nanmax(age_all)))
    (vdir / "index.html").write_text(
        build_html(results, age_range, n_members, per_method_pngs),
        encoding="utf-8")

    print("[validate] metrics:")
    for k, v in results.items():
        print(f"  {k}: r={v.get('r'):.3f} bias={v.get('bias'):+.3f} maxD={v.get('maxD'):.3f} midHol={v.get('midHol_mine'):.2f}|{v.get('midHol_pub'):.2f}")
    print(f"[validate] wrote {vdir}/index.html + {len(per_method_pngs)} per-method PNGs + combined overlay")


if __name__ == "__main__":
    main()
