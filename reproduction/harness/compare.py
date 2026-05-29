#!/usr/bin/env python3
"""Compare a harness per-method CSV (_repro/out/<m>_global.csv) to published,
anchored to 0 at 100 BP. Prints the same metrics as _audit/fidelity.py plus the
prior BAM/nens=50 snapshot for that method, so the reproduction gain is visible.
"""
import argparse
import numpy as np
import pandas as pd

PUB = "presto-temp12k/reference_data/published"
TABLE1 = {"scc": 0.50, "dcc": 0.50, "gam": 0.44, "cps": 1.08, "paico": 0.42}
# prior snapshot (BAM ages, pkl, template filter, nens=50) for reference
PRIOR = {"dcc": dict(r=0.998, bias=0.043, rmse=0.053, maxd=0.092, amp=1.07, mh=0.56, c12=-0.80, spread=0.97),
         "scc": dict(r=0.991, bias=0.061, rmse=0.076, maxd=0.161, amp=1.06, mh=0.60, c12=-0.75, spread=0.79),
         "gam": dict(r=0.969, bias=0.085, rmse=0.105, maxd=0.260, amp=1.11, mh=0.53, c12=-0.59, spread=0.22),
         "paico": dict(r=0.998, bias=0.037, rmse=0.091, maxd=0.294, amp=1.29, mh=0.54, c12=-0.95, spread=1.08),
         "cps": dict(r=0.996, bias=0.235, rmse=0.255, maxd=0.525, amp=0.96, mh=1.29, c12=-2.97, spread=1.00)}


def anchor(age, y, at=100.0):
    o = np.argsort(age)
    return y - np.interp(at, age[o], y[o])


def wmean(age, y, lo, hi):
    m = (age >= lo) & (age <= hi)
    return float(np.nanmean(y[m])) if m.any() else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True)
    ap.add_argument("--csv", required=True)
    a = ap.parse_args()
    m = a.method

    d = pd.read_csv(a.csv)
    age = d["binAges"].to_numpy(float)
    ens = d.drop(columns=["binAges"]).to_numpy(float)
    o = np.argsort(age); age = age[o]; ens = ens[o]
    med = np.nanmedian(ens, 1); q05 = np.nanpercentile(ens, 5, 1); q95 = np.nanpercentile(ens, 95, 1)

    p = pd.read_csv(f"{PUB}/{m}_published.csv")
    pa = p["age_bp"].to_numpy(float); po = np.argsort(pa)
    pmed = np.interp(age, pa[po], p["median"].to_numpy(float)[po])
    plo = np.interp(age, pa[po], p["lo"].to_numpy(float)[po])
    phi = np.interp(age, pa[po], p["hi"].to_numpy(float)[po])

    me = anchor(age, med); pu = anchor(age, pmed)
    msk = np.isfinite(me) & np.isfinite(pu)
    r = float(np.corrcoef(me[msk], pu[msk])[0, 1])
    bias = float(np.nanmean(me - pu)); rmse = float(np.sqrt(np.nanmean((me - pu) ** 2)))
    maxd = float(np.nanmax(np.abs(me - pu))); amp = float(np.nanstd(me) / np.nanstd(pu))
    mh = wmean(age, me, 5500, 6500); mhp = wmean(age, pu, 5500, 6500)
    c12 = wmean(age, me, 11500, 12000); c12p = wmean(age, pu, 11500, 12000)
    spread = float(np.nanmean(q95 - q05) / np.nanmean(phi - plo))
    flat = max((1 + sum(1 for j in range(i, len(med)) if med[j] == med[i]) for i in []), default=1)
    # simple flat-run
    best = run = 1
    for i in range(1, len(med)):
        run = run + 1 if med[i] == med[i - 1] else 1
        best = max(best, run)

    pr = PRIOR.get(m, {})
    print(f"\n=== {m.upper()} reproduction (real ageEnsemble + temp12kEnsemble) vs prior (BAM/pkl/template) ===")
    print(f"{'metric':<10}{'NOW':>10}{'PRIOR':>10}{'published/target':>20}")
    rows = [("r", r, pr.get("r"), "1.000"),
            ("bias", bias, pr.get("bias"), "0.000"),
            ("RMSE", rmse, pr.get("rmse"), "0.000"),
            ("maxD", maxd, pr.get("maxd"), "0.000"),
            ("amp", amp, pr.get("amp"), "1.00"),
            ("midHol", mh, pr.get("mh"), f"{mhp:.2f} (T1 {TABLE1.get(m,'?')})"),
            ("12ka", c12, pr.get("c12"), f"{c12p:.2f}"),
            ("spread", spread, pr.get("spread"), "1.00"),
            ("flat", best, None, "1")]
    for name, now, prior, tgt in rows:
        ps = f"{prior:.3f}" if isinstance(prior, (int, float)) else "-"
        nv = f"{now:.3f}" if isinstance(now, float) else str(now)
        print(f"{name:<10}{nv:>10}{ps:>10}{tgt:>20}")


if __name__ == "__main__":
    main()
