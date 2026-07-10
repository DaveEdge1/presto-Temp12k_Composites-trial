#!/usr/bin/env python3
"""Standalone re-impl of reproduction/harness/compare.py (no pandas).
Anchors both curves to 0 at 100 BP, prints the same metrics."""
import argparse, csv, numpy as np

import os
PUB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "reference_data", "published")
TABLE1 = {"scc":0.50,"dcc":0.50,"gam":0.44,"cps":1.08,"paico":0.42}

def load_csv(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    header = rows[0]
    data = np.array([[float(x) if x not in ("","NA") else np.nan for x in r] for r in rows[1:]], float)
    return header, data

def anchor(age, y, at=100.0):
    o = np.argsort(age); return y - np.interp(at, age[o], y[o])

def wmean(age, y, lo, hi):
    m = (age>=lo)&(age<=hi); return float(np.nanmean(y[m])) if m.any() else float("nan")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--method",required=True); ap.add_argument("--csv",required=True)
    a = ap.parse_args(); m = a.method
    h, d = load_csv(a.csv)
    ai = h.index("binAges"); age = d[:,ai]
    ens = np.delete(d, ai, axis=1)
    o = np.argsort(age); age = age[o]; ens = ens[o]
    med = np.nanmedian(ens,1); q05 = np.nanpercentile(ens,5,1); q95 = np.nanpercentile(ens,95,1)
    ph, pd_ = load_csv(f"{PUB}/{m}_published.csv")
    pa = pd_[:,ph.index("age_bp")]; po = np.argsort(pa)
    pmed = np.interp(age, pa[po], pd_[:,ph.index("median")][po])
    plo  = np.interp(age, pa[po], pd_[:,ph.index("lo")][po])
    phi  = np.interp(age, pa[po], pd_[:,ph.index("hi")][po])
    me = anchor(age, med); pu = anchor(age, pmed)
    msk = np.isfinite(me)&np.isfinite(pu)
    r = float(np.corrcoef(me[msk],pu[msk])[0,1])
    bias = float(np.nanmean(me-pu)); rmse = float(np.sqrt(np.nanmean((me-pu)**2)))
    maxd = float(np.nanmax(np.abs(me-pu))); amp = float(np.nanstd(me)/np.nanstd(pu))
    mh = wmean(age,me,5500,6500); mhp = wmean(age,pu,5500,6500)
    c12 = wmean(age,me,11500,12000); c12p = wmean(age,pu,11500,12000)
    spread = float(np.nanmean(q95-q05)/np.nanmean(phi-plo))
    print(f"\n=== {m.upper()} harness vs published (anchored 100 BP) ===")
    for name,val,tgt in [("r",r,"1.000"),("bias",bias,"0.000"),("RMSE",rmse,"0.000"),
        ("maxD",maxd,"0.000"),("amp",amp,"1.00"),("midHol",mh,f"{mhp:.2f} (T1 {TABLE1.get(m,'?')})"),
        ("12ka",c12,f"{c12p:.2f}"),("spread",spread,"1.00")]:
        print(f"{name:<8}{val:>10.3f}{tgt:>22}")

if __name__ == "__main__": main()
