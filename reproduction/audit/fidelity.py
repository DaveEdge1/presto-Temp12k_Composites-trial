#!/usr/bin/env python3
"""Fidelity audit: my reconstruction vs published Kaufman 2020 curves.

Both curves are put on the OFFICIAL footing (readme): each method's median
anchored to 0 at 100 yr BP. Then absolute (not recentered) comparison so that
amplitude and offset errors are visible — unlike the r/CE-on-recentered metric,
which hides them.
"""
import glob
import numpy as np
import pandas as pd

RECON = "_audit/recon"
PUB = "presto-temp12k/reference_data/published"
METHODS = ["scc", "dcc", "gam", "cps", "paico"]
TABLE1 = {"scc": 0.50, "dcc": 0.50, "gam": 0.44, "cps": 1.08, "paico": 0.42,
          "consensus": 0.70}  # mid-Holocene 6.5-5.5ka, paper Table 1 / abstract


def anchor(age, y, at=100.0):
    """Subtract the value at `at` BP (linear interp) so y(at)=0."""
    return y - np.interp(at, age[np.argsort(age)], y[np.argsort(age)])


def my_method_curve(m):
    """median, q05, q95 from my per-method global ensemble."""
    df = pd.read_csv(f"{RECON}/methods/{m}_global.csv")
    age = df["binAges"].to_numpy(float)
    ens = df.drop(columns=["binAges"]).to_numpy(float)
    o = np.argsort(age)
    return (age[o], np.nanmedian(ens, axis=1)[o],
            np.nanpercentile(ens, 5, axis=1)[o],
            np.nanpercentile(ens, 95, axis=1)[o])


def my_consensus_curve():
    df = pd.read_csv(f"{RECON}/reconstruction.csv")
    age = df["age_bp"].to_numpy(float)
    o = np.argsort(age)
    return (age[o], df["mean"].to_numpy(float)[o],
            df["lo_95"].to_numpy(float)[o], df["hi_95"].to_numpy(float)[o])


def pub_curve(m):
    name = "consensus" if m == "consensus" else m
    df = pd.read_csv(f"{PUB}/{name}_published.csv")
    age = df["age_bp"].to_numpy(float)
    o = np.argsort(age)
    return (age[o], df["median"].to_numpy(float)[o],
            df["lo"].to_numpy(float)[o], df["hi"].to_numpy(float)[o])


def window_mean(age, y, lo, hi):
    m = (age >= lo) & (age <= hi)
    return float(np.nanmean(y[m])) if m.any() else float("nan")


def clip_run(y):
    """Longest run of identical consecutive values (flat-clip detector)."""
    best = run = 1
    for i in range(1, len(y)):
        run = run + 1 if y[i] == y[i - 1] else 1
        best = max(best, run)
    return best


print(f"{'method':<10}{'r':>7}{'bias':>8}{'RMSE':>8}{'maxD':>8}"
      f"{'amp(mine/pub)':>15}{'midHol m|p|T1':>20}{'12ka m|p':>14}"
      f"{'spread m/p':>12}{'flat':>6}")
print("-" * 118)

rows = []
for m in METHODS + ["consensus"]:
    if m == "consensus":
        a, med, q05, q95 = my_consensus_curve()
    else:
        a, med, q05, q95 = my_method_curve(m)
    pa, pmed, plo, phi = pub_curve(m)

    # common grid = my ages (published is on same 100-yr grid)
    pmed_i = np.interp(a, pa, pmed)
    plo_i = np.interp(a, pa, plo)
    phi_i = np.interp(a, pa, phi)

    # anchor both to 0 at 100 BP (official convention)
    me = anchor(a, med)
    pu = anchor(a, pmed_i)

    msk = np.isfinite(me) & np.isfinite(pu)
    r = float(np.corrcoef(me[msk], pu[msk])[0, 1])
    bias = float(np.nanmean(me - pu))
    rmse = float(np.sqrt(np.nanmean((me - pu) ** 2)))
    maxd = float(np.nanmax(np.abs(me - pu)))
    amp = float(np.nanstd(me) / np.nanstd(pu))
    mh_m = window_mean(a, me, 5500, 6500)
    mh_p = window_mean(a, pu, 5500, 6500)
    c12_m = window_mean(a, me, 11500, 12000)
    c12_p = window_mean(a, pu, 11500, 12000)
    spread_m = float(np.nanmean(q95 - q05))
    spread_p = float(np.nanmean(phi_i - plo_i))
    spread_ratio = spread_m / spread_p if spread_p else float("nan")
    flat = clip_run(med)

    print(f"{m:<10}{r:>7.3f}{bias:>8.3f}{rmse:>8.3f}{maxd:>8.3f}"
          f"{amp:>15.2f}{mh_m:>7.2f}|{mh_p:.2f}|{TABLE1[m]:.2f}"
          f"{c12_m:>8.2f}|{c12_p:.2f}{spread_ratio:>12.2f}{flat:>6d}")
    rows.append(m)

print("\nNotes:")
print("  amp  = std(mine)/std(pub) over full record (1.0 = matched amplitude)")
print("  midHol = mean over 5.5-6.5 ka (mine | published | paper Table 1)")
print("  12ka = mean over 11.5-12 ka (mine | published) — deglacial cold end")
print("  spread = mean(q95-q05): mine/published ratio (1.0 = matched uncertainty)")
print("  flat = longest run of identical consecutive median values (clipping flag)")
