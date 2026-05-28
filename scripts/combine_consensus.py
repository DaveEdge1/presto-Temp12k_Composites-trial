#!/usr/bin/env python3
"""Combine per-method global ensembles into the multi-method consensus GMST.

Reads each enabled method's results/methods/<m>_global.csv (binAges + ensemble
columns), pools ALL members across methods into a single multi-method ensemble
(the paper's 2500-member consensus = 5 methods x 500), and writes:

  results/reconstruction.csv
      year, age_bp, mean (=consensus median), lo_95 (5th), hi_95 (95th),
      n_members, and one <method>_median column per contributing method.

This is the headline product (Kaufman et al. 2020, Fig. 3). The per-method
medians are kept for the overlay figure. Downstream outputs_to_netcdf.py +
make_figures.py consume reconstruction.csv.
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods-dir", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.methods_dir, "*_global.csv")))
    if not files:
        raise SystemExit(f"[combine] no *_global.csv in {args.methods_dir}")

    bin_ages = None
    pooled = []          # list of (nbins x nmembers) arrays
    per_method_median = {}
    for f in files:
        method = os.path.basename(f).split("_")[0]
        df = pd.read_csv(f)
        ba = df["binAges"].to_numpy(dtype=float)
        ens = df.drop(columns=["binAges"]).to_numpy(dtype=float)
        if bin_ages is None:
            bin_ages = ba
        elif not np.allclose(ba, bin_ages, equal_nan=True):
            # reindex onto the first method's binAges if grids differ
            idx = [int(np.argmin(np.abs(bin_ages - a))) for a in ba]
            tmp = np.full((bin_ages.size, ens.shape[1]), np.nan)
            tmp[idx, :] = ens
            ens = tmp
        pooled.append(ens)
        per_method_median[method] = np.nanmedian(ens, axis=1)
        print(f"[combine] {method}: {ens.shape[1]} members")

    allens = np.hstack(pooled)                          # nbins x total_members
    finite_per_bin = np.sum(np.isfinite(allens), axis=1)
    with np.errstate(all="ignore"):
        med = np.nanmedian(allens, axis=1)
        lo = np.nanpercentile(allens, 5, axis=1)
        hi = np.nanpercentile(allens, 95, axis=1)

    out = pd.DataFrame({
        "year": (1950 - bin_ages).astype(int),          # CE (negative = BCE)
        "age_bp": bin_ages,
        "mean": med,                                     # consensus = ensemble median
        "lo_95": lo,
        "hi_95": hi,
        "n_members": finite_per_bin,
    })
    for m, v in per_method_median.items():
        out[f"{m}_median"] = v

    out = out.sort_values("year").reset_index(drop=True)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"[combine] wrote {args.out_csv}: {len(out)} bins, "
          f"{allens.shape[1]} pooled members across {len(files)} methods")


if __name__ == "__main__":
    main()
