#!/usr/bin/env python3
"""GAM stage 2: read pooled (age, temp, cell, band) samples, fit one LinearGAM
per cell on the pooled ensemble cloud (matches the published gam_ensemble.py),
draw n_draws posterior samples per cell, aggregate cells->band (equal-area grid
mean) and bands->global (sin-area weights), apply the archival reference.

Output: gam_global.csv (binAges + nens columns) matching the R harness schema.
"""
import argparse
import sys
import numpy as np
import pandas as pd
from pygam import LinearGAM, s

ap = argparse.ArgumentParser()
ap.add_argument("--inp", default="/repro/out/gam_pooled.csv")
ap.add_argument("--out", default="/repro/out/gam_global.csv")
ap.add_argument("--nens", type=int, default=50)
ap.add_argument("--ncores", type=int, default=12)
a = ap.parse_args()

LATBINS = np.arange(-90, 91, 30)
N_BANDS = 6
ZW = np.sin(LATBINS[1:] * np.pi / 180) - np.sin(LATBINS[:-1] * np.pi / 180); ZW /= ZW.sum()
binvec = np.arange(-50, 12051, 100, dtype=float)
bin_ages = (binvec[:-1] + binvec[1:]) / 2.0
NB = bin_ages.size

print(f"[gam] loading {a.inp} ...", flush=True)
d = pd.read_csv(a.inp)
print(f"[gam] {len(d):,} points, {d['cell'].nunique()} cells, "
      f"age {d['age'].min():.0f}-{d['age'].max():.0f}", flush=True)


def fit_cell(grp_xy_band):
    cid, x, y, band = grp_xy_band
    if x.size < 20 or np.ptp(x) < 200:
        return cid, None, band
    # per-cell anchor to 3-5 ka (paper's GAM modern_old=5000, modern_young=3000)
    ref = (x >= 3000) & (x <= 5000)
    y0 = y - (np.nanmean(y[ref]) if ref.any() else np.nanmean(y))
    try:
        gam = LinearGAM(s(0)).gridsearch(x[:, None], y0, progress=False)
        # posterior draws of the spline coefficients evaluated at bin_ages
        draws = gam.sample(x[:, None], y0, sample_at_X=bin_ages[:, None],
                           n_draws=a.nens, n_bootstraps=1)
        # mask outside the cell's age coverage so we don't extrapolate
        outside = (bin_ages < x.min()) | (bin_ages > x.max())
        draws[:, outside] = np.nan
        return cid, draws.T, band                     # (NB, nens)
    except Exception as exc:
        print(f"[gam] cell {cid}: {exc}", file=sys.stderr)
        return cid, None, band


# Pre-extract per-cell numpy arrays, then fit in parallel
groups = []
for cid, grp in d.groupby("cell"):
    groups.append((int(cid), grp["age"].to_numpy(float),
                   grp["temp"].to_numpy(float), int(grp["band"].iloc[0])))
print(f"[gam] fitting {len(groups)} cells on {a.ncores} cores ...", flush=True)

import multiprocessing as mp
with mp.get_context("fork").Pool(a.ncores) as pool:
    results = pool.map(fit_cell, groups)

cell_draws = {cid: dr for cid, dr, _ in results if dr is not None}
cell_band = {cid: bn for cid, _, bn in results}
print(f"[gam] fit {len(cell_draws)} cells", flush=True)

# Per band: simple mean of cells in band (equal-area grid -> no extra weights)
band_ens = np.full((NB, N_BANDS, a.nens), np.nan)
for b in range(N_BANDS):
    cells_b = [c for c, bn in cell_band.items()
               if bn == b + 1 and c in cell_draws]
    if not cells_b:
        continue
    arr = np.stack([cell_draws[c] for c in cells_b], axis=0)   # (n_cells, NB, nens)
    band_ens[:, b, :] = np.nanmean(arr, axis=0)

# Per draw: area-weight bands -> global
glob = np.full((NB, a.nens), np.nan)
for k in range(a.nens):
    bm = band_ens[:, :, k]
    w = np.tile(ZW, (NB, 1)).astype(float); w[~np.isfinite(bm)] = np.nan
    num = np.nansum(bm * w, axis=1); den = np.nansum(w, axis=1)
    glob[:, k] = np.where(den > 0, num / den, np.nan)

# Archival reference: per-member subtract full-12k mean, then median=0 at 100 BP
glob = glob - np.nanmean(glob, axis=0, keepdims=True)
r100 = int(np.argmin(np.abs(bin_ages - 100)))
glob = glob - np.nanmedian(glob[r100, :])

df = pd.DataFrame(glob, columns=[f"ens{i+1}" for i in range(a.nens)])
df.insert(0, "binAges", bin_ages)
df.to_csv(a.out, index=False)
print(f"[gam] wrote {a.out}", flush=True)
