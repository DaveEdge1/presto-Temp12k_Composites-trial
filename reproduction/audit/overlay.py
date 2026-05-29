#!/usr/bin/env python3
"""Per-method overlay: mine (median + 5-95) vs published, anchored to 0 at 100 BP."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RECON = "_audit/recon"
PUB = "presto-temp12k/reference_data/published"
METHODS = ["scc", "dcc", "gam", "cps", "paico", "consensus"]


def anchor(age, y, at=100.0):
    o = np.argsort(age)
    return y - np.interp(at, age[o], y[o])


def mine(m):
    if m == "consensus":
        d = pd.read_csv(f"{RECON}/reconstruction.csv")
        a = d["age_bp"].to_numpy(float)
        return a, d["mean"].to_numpy(float), d["lo_95"].to_numpy(float), d["hi_95"].to_numpy(float)
    d = pd.read_csv(f"{RECON}/methods/{m}_global.csv")
    a = d["binAges"].to_numpy(float)
    ens = d.drop(columns=["binAges"]).to_numpy(float)
    return a, np.nanmedian(ens, 1), np.nanpercentile(ens, 5, 1), np.nanpercentile(ens, 95, 1)


def pub(m):
    d = pd.read_csv(f"{PUB}/{m}_published.csv")
    return (d["age_bp"].to_numpy(float), d["median"].to_numpy(float),
            d["lo"].to_numpy(float), d["hi"].to_numpy(float))


fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
for ax, m in zip(axes.ravel(), METHODS):
    a, md, lo, hi = mine(m)
    pa, pmd, plo, phi = pub(m)
    a0 = a / 1000.0; pa0 = pa / 1000.0
    off_m = np.interp(100.0, a[np.argsort(a)], md[np.argsort(a)])
    off_p = np.interp(100.0, pa[np.argsort(pa)], pmd[np.argsort(pa)])
    ax.fill_between(a0, lo - off_m, hi - off_m, color="tab:blue", alpha=0.2, label="mine 5-95")
    ax.plot(a0, md - off_m, color="tab:blue", lw=2, label="mine median")
    ax.fill_between(pa0, plo - off_p, phi - off_p, color="tab:red", alpha=0.15, label="pub 5-95")
    ax.plot(pa0, pmd - off_p, color="tab:red", lw=1.6, ls="--", label="published median")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlim(12, 0); ax.set_title(m.upper()); ax.grid(alpha=0.3, lw=0.5)
    ax.set_xlabel("ka BP"); ax.set_ylabel("anomaly (degC)")
axes.ravel()[0].legend(fontsize=7, loc="lower center")
fig.suptitle("Temp12k fidelity: reconstruction (blue) vs published Kaufman 2020 (red), anchored at 100 BP", y=1.0)
fig.tight_layout()
fig.savefig("_audit/fidelity_overlay.png", dpi=130, bbox_inches="tight")
print("wrote _audit/fidelity_overlay.png")
