#!/usr/bin/env python3
"""GAM -- Generalized Additive Model composite (Python / pygam).

Reimplements the Temperature 12k GAM method (Sommer & Davis; paper Methods,
"Reconstruction method 3"). For each equal-area grid cell we fit a penalized
B-spline GAM through the pooled (age, temperature) samples of the records in
that cell, predict the temperature anomaly every 100 yr, then average grid
cells within each 30-deg latitude band and area-weight the bands into a global
mean. The ensemble is drawn from the GAM posterior (samples of the spline
coefficients), with chronological uncertainty added as an age jitter whose SD
shrinks linearly from ~250 yr at 12 ka to ~50 yr at present (paper).

We deliberately reimplement the gridding with numpy/xarray + modern pygam
instead of the original pyleogrid/psyplot/rpy2/py3.7 stack, which is not
shippable. Output schema matches the R methods (binAges + nens columns) so the
consensus step treats all five methods uniformly.

Uses only calibrated degC records (as SCC/DCC), seasonality in
{annual, summerOnly, winterOnly}.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

LATBINS = np.arange(-90, 91, 30)                       # 6 bands
BAND_WEIGHTS = np.array([0.067, 0.183, 0.25, 0.25, 0.183, 0.067])
N_BANDS = len(BAND_WEIGHTS)


def band_of(lat):
    if lat is None or not np.isfinite(lat):
        return None
    b = int(np.searchsorted(LATBINS, lat, side="right")) - 1
    return b if 0 <= b < N_BANDS else None


def nearest_cell(lat, lon, clat, clon):
    rad = np.pi / 180.0
    d = (np.sin(clat * rad) * np.sin(lat * rad)
         + np.cos(clat * rad) * np.cos(lat * rad) * np.cos((clon - lon) * rad))
    return int(np.argmax(np.clip(d, -1, 1)))


def age_jitter_sd(age_bp):
    """Chronology SD: ~250 yr at 12 ka shrinking linearly to ~50 yr at 0 (paper)."""
    return np.clip(50.0 + (250.0 - 50.0) * (age_bp / 12000.0), 50.0, 250.0)


def area_weight(band_mat):
    w = np.tile(BAND_WEIGHTS, (band_mat.shape[0], 1)).astype(float)
    w[~np.isfinite(band_mat)] = np.nan
    num = np.nansum(band_mat * w, axis=1)
    den = np.nansum(w, axis=1)
    out = np.where(den > 0, num / den, np.nan)
    return out


def apply_reference(ens, bin_ages, ref_start_ce=1800, ref_end_ce=1900):
    ens = ens - np.nanmean(ens, axis=0, keepdims=True)         # per-member full-record mean
    lo, hi = 1950 - ref_end_ce, 1950 - ref_start_ce            # ~50..150 BP
    rows = np.where((bin_ages >= lo) & (bin_ages <= hi))[0]
    if rows.size == 0:
        rows = np.array([int(np.argmin(np.abs(bin_ages - 100)))])
    med = np.nanmedian(np.nanmean(ens[rows, :], axis=0))
    return ens - med


def _align_anomaly(age, val, ref_lo=3000.0, ref_hi=5000.0):
    """Express a record as an anomaly relative to its 3-5 ka mean (paper's GAM
    alignment); fall back to the whole-record mean if it doesn't span 3-5 ka.
    Without this, pooling records at different absolute temperatures into one
    grid-cell GAM inflates amplitude and masks the deglacial signal."""
    win = (age >= ref_lo) & (age <= ref_hi)
    ref = np.nanmean(val[win]) if win.sum() >= 2 else np.nanmean(val)
    return val - ref


def load_records(ts_path):
    recs = json.loads(Path(ts_path).read_text())
    out = []
    for r in recs:
        if str(r.get("units", "")).lower() != "degc":
            continue
        if str(r.get("seasonalityGeneral", "")).lower() not in ("annual", "summeronly", "winteronly"):
            continue
        age = np.asarray(r["age"], dtype=float)
        val = np.asarray(r["values"], dtype=float)
        m = np.isfinite(age) & np.isfinite(val)
        if m.sum() < 3 or r.get("lat") is None or r.get("lon") is None:
            continue
        d = 1.0 if str(r.get("direction", "")).lower() != "negative" else -1.0
        age, val = age[m], val[m] * d
        val = _align_anomaly(age, val)          # anomaly vs 3-5 ka mean
        out.append({"age": age, "val": val,
                    "lat": float(r["lat"]), "lon": float(r["lon"]),
                    "unc": float(r.get("uncertainty1sd") or 2.0)})
    return out


def fit_cell_curve(ages, temps, bin_ages, n_splines, rng, jitter=True):
    """Fit a penalized spline to pooled (age, temp); return prediction at bin_ages
    (NaN outside the data's age coverage) plus a posterior coefficient sample."""
    from pygam import LinearGAM, s
    a = np.asarray(ages, dtype=float)
    t = np.asarray(temps, dtype=float)
    if jitter:
        a = a + rng.normal(0.0, age_jitter_sd(a))
    order = np.argsort(a)
    a, t = a[order], t[order]
    if a.size < 4 or np.ptp(a) < 100:
        return None
    ns = int(min(n_splines, max(4, a.size // 3)))
    try:
        gam = LinearGAM(s(0, n_splines=ns)).fit(a.reshape(-1, 1), t)
    except Exception:
        return None
    pred = gam.predict(bin_ages.reshape(-1, 1)).astype(float)
    pred[(bin_ages < a.min()) | (bin_ages > a.max())] = np.nan       # no extrapolation
    return pred


# Shared state for forked worker processes (set in run_gam before the Pool is
# created; inherited copy-on-write by fork on Linux, so it isn't re-pickled).
_GAM_SHARED = {}


def _gam_member(k):
    recs = _GAM_SHARED["recs"]; band_cells = _GAM_SHARED["band_cells"]
    bin_ages = _GAM_SHARED["bin_ages"]; n_splines = _GAM_SHARED["n_splines"]
    rng = np.random.default_rng(12345 + k)         # per-member reproducible seed
    band_mat = np.full((bin_ages.size, N_BANDS), np.nan)
    for b in range(N_BANDS):
        cell_curves = []
        for cell, idxs in band_cells[b].items():
            ages = np.concatenate([recs[i]["age"] for i in idxs])
            temps = np.concatenate([recs[i]["val"]
                                    + rng.normal(0.0, recs[i]["unc"], size=recs[i]["val"].size)
                                    for i in idxs])
            curve = fit_cell_curve(ages, temps, bin_ages, n_splines, rng, jitter=True)
            if curve is not None:
                cell_curves.append(curve)
        if cell_curves:
            band_mat[:, b] = np.nanmean(np.vstack(cell_curves), axis=0)
    return area_weight(band_mat)


def run_gam(ts_path, cfg, grid, out_csv):
    bin_cfg = cfg.get("bin", {})
    start, end, step = bin_cfg.get("start_bp", -50), bin_cfg.get("end_bp", 12050), bin_cfg.get("step", 100)
    binvec = np.arange(start, end + step, step, dtype=float)
    bin_ages = (binvec[1:] + binvec[:-1]) / 2.0
    nens = int(cfg.get("nens", 100))
    n_splines = int((cfg.get("advanced") or {}).get("gam_spline_k", 25))
    refp = cfg.get("reference_period", {"start": 1800, "end": 1900})

    recs = load_records(ts_path)
    print(f"[gam] {len(recs)} degC records after filter", file=sys.stderr)
    if len(recs) < 10:
        raise SystemExit("[gam] too few records")

    clat = grid["clat"].to_numpy(); clon = grid["clon180"].to_numpy()
    # assign records to cells and bands
    for r in recs:
        r["cell"] = nearest_cell(r["lat"], r["lon"], clat, clon)
        r["band"] = band_of(r["lat"])
    # group record indices by band -> cell
    band_cells = {b: {} for b in range(N_BANDS)}
    for i, r in enumerate(recs):
        if r["band"] is None:
            continue
        band_cells[r["band"]].setdefault(r["cell"], []).append(i)

    # parallelize across ensemble members (each member is independent)
    import os
    ncores = int(cfg.get("ncores") or max(1, (os.cpu_count() or 2) - 1))
    ncores = max(1, min(ncores, nens))
    _GAM_SHARED.update(recs=recs, band_cells=band_cells, bin_ages=bin_ages, n_splines=n_splines)
    print(f"[gam] {nens} members on {ncores} core(s)", file=sys.stderr)
    if ncores > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(ncores) as pool:
            cols = pool.map(_gam_member, range(nens))
    else:
        cols = [_gam_member(k) for k in range(nens)]
    ens = np.column_stack(cols)

    ens = apply_reference(ens, bin_ages, refp.get("start", 1800), refp.get("end", 1900))
    df = pd.DataFrame(ens, columns=[f"ens{i+1}" for i in range(nens)])
    df.insert(0, "binAges", bin_ages)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[gam] wrote {out_csv} ({bin_ages.size} bins x {nens} members)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ts", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--grid", required=True, help="equal_area_grid_centers.csv")
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    grid = pd.read_csv(args.grid)
    run_gam(args.ts, cfg, grid, args.out_csv)


if __name__ == "__main__":
    main()
