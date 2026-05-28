#!/usr/bin/env python3
"""Pool the per-method ensembles into the multi-method consensus.

Reads results/methods/<m>_global.csv (binAges + ensemble columns) and
<m>_bands.csv (long: binAges, band, ensemble columns), pools ALL members
across methods (the paper's multi-method ensemble = 5 methods x nens), and writes:

  reconstruction.csv          year, age_bp, mean(=median), lo_95, hi_95, n_members,
                              and one <method>_median column per method (headline GMST)
  reconstruction_ensemble.csv binAges + pooled global ensemble members (-> tas_gm)
  reconstruction_zonal.csv    binAges + per-30deg-band consensus (median/lo/hi) (-> tas stripes)

Downstream: outputs_to_netcdf.py builds the spatial NetCDF, make_figures.py / validate.py.
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

N_BANDS = 6


def _reindex(ens, ba, bin_ages):
    if np.allclose(ba, bin_ages, equal_nan=True):
        return ens
    idx = [int(np.argmin(np.abs(bin_ages - a))) for a in ba]
    tmp = np.full((bin_ages.size, ens.shape[1]), np.nan)
    tmp[idx, :] = ens
    return tmp


def pool_global(methods_dir):
    files = sorted(glob.glob(os.path.join(methods_dir, "*_global.csv")))
    if not files:
        raise SystemExit(f"[combine] no *_global.csv in {methods_dir}")
    bin_ages, pooled, per_method = None, [], {}
    for f in files:
        method = os.path.basename(f).split("_")[0]
        df = pd.read_csv(f)
        ba = df["binAges"].to_numpy(float)
        ens = df.drop(columns=["binAges"]).to_numpy(float)
        if bin_ages is None:
            bin_ages = ba
        ens = _reindex(ens, ba, bin_ages)
        pooled.append(ens)
        per_method[method] = np.nanmedian(ens, axis=1)
        print(f"[combine] {method}: {ens.shape[1]} global members")
    return bin_ages, np.hstack(pooled), per_method


def pool_zonal(methods_dir, bin_ages):
    files = sorted(glob.glob(os.path.join(methods_dir, "*_bands.csv")))
    if not files:
        return None
    band_pool = {b: [] for b in range(1, N_BANDS + 1)}
    for f in files:
        df = pd.read_csv(f)
        for b in range(1, N_BANDS + 1):
            sub = df[df["band"] == b].sort_values("binAges")
            if sub.empty:
                continue
            ba = sub["binAges"].to_numpy(float)
            ens = sub.drop(columns=["binAges", "band"]).to_numpy(float)
            band_pool[b].append(_reindex(ens, ba, bin_ages))
    zonal = {"binAges": bin_ages}
    with np.errstate(all="ignore"):
        for b in range(1, N_BANDS + 1):
            if band_pool[b]:
                allb = np.hstack(band_pool[b])
                zonal[f"band{b}_med"] = np.nanmedian(allb, axis=1)
                zonal[f"band{b}_lo"] = np.nanpercentile(allb, 5, axis=1)
                zonal[f"band{b}_hi"] = np.nanpercentile(allb, 95, axis=1)
            else:
                nan = np.full(bin_ages.size, np.nan)
                zonal[f"band{b}_med"] = zonal[f"band{b}_lo"] = zonal[f"band{b}_hi"] = nan
    return pd.DataFrame(zonal)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods-dir", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()
    out_dir = Path(args.out_csv).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    bin_ages, allens, per_method = pool_global(args.methods_dir)
    n_members = np.sum(np.isfinite(allens), axis=1)
    with np.errstate(all="ignore"):
        med = np.nanmedian(allens, axis=1)
        lo = np.nanpercentile(allens, 5, axis=1)
        hi = np.nanpercentile(allens, 95, axis=1)

    out = pd.DataFrame({
        "year": (1950 - bin_ages).astype(int),
        "age_bp": bin_ages, "mean": med, "lo_95": lo, "hi_95": hi,
        "n_members": n_members,
    })
    for m, v in per_method.items():
        out[f"{m}_median"] = v
    out = out.sort_values("year").reset_index(drop=True)
    out.to_csv(args.out_csv, index=False)
    print(f"[combine] wrote {args.out_csv}: {len(out)} bins, "
          f"{allens.shape[1]} pooled members across {len(per_method)} methods")

    # pooled global ensemble (for tas_gm)
    ens_df = pd.DataFrame(allens, columns=[f"ens{i+1}" for i in range(allens.shape[1])])
    ens_df.insert(0, "binAges", bin_ages)
    ens_path = out_dir / "reconstruction_ensemble.csv"
    ens_df.to_csv(ens_path, index=False)
    print(f"[combine] wrote {ens_path} ({allens.shape[1]} pooled global members)")

    # zonal consensus per band (for tas stripes)
    zonal = pool_zonal(args.methods_dir, bin_ages)
    if zonal is not None:
        zonal_path = out_dir / "reconstruction_zonal.csv"
        zonal.to_csv(zonal_path, index=False)
        print(f"[combine] wrote {zonal_path} (per-band consensus, {N_BANDS} bands)")


if __name__ == "__main__":
    main()
