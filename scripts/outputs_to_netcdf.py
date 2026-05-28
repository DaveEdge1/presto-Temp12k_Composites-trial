#!/usr/bin/env python3
"""reconstruction.csv -> 1D CF-NetCDF (consensus GMST + per-method medians).

1D output (time only) -> visualize.yml routes to the static-Pages tile UI.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def convert(in_csv: Path, out_nc: Path) -> None:
    df = pd.read_csv(in_csv).sort_values("year").reset_index(drop=True)

    def col(name):
        return df[name].to_numpy(dtype=float) if name in df.columns \
            else np.full(len(df), np.nan)

    data = {
        "gmst":        ("time", col("mean")),
        "gmst_lo_95":  ("time", col("lo_95")),
        "gmst_hi_95":  ("time", col("hi_95")),
        "n_members":   ("time", col("n_members")),
    }
    for m in [c for c in df.columns if c.endswith("_median")]:
        data[m] = ("time", col(m))

    ds = xr.Dataset(
        data_vars=data,
        coords={"time": df["year"].to_numpy(dtype=int)},
        attrs={
            "title": "Temperature 12k multi-method composite GMST",
            "source": "presto-temp12k (SCC/DCC/GAM/CPS/PaiCo consensus)",
            "reference": "Kaufman et al. 2020, Sci. Data 7:201, doi:10.1038/s41597-020-0530-7",
            "Conventions": "CF-1.10",
        },
    )
    ds["time"].attrs.update(units="years_AD", standard_name="time", long_name="year CE")
    ds["gmst"].attrs.update(long_name="consensus GMST anomaly", units="degC")
    ds["gmst_lo_95"].attrs.update(long_name="5th percentile", units="degC")
    ds["gmst_hi_95"].attrs.update(long_name="95th percentile", units="degC")

    out_nc.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_nc)
    print(f"[outputs_to_netcdf] wrote {out_nc} ({len(df)} time steps)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", required=True, type=Path)
    ap.add_argument("--out-nc", required=True, type=Path)
    args = ap.parse_args()
    convert(args.in_csv, args.out_nc)


if __name__ == "__main__":
    main()
