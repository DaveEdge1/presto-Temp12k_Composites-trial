#!/usr/bin/env python3
"""reconstruction.csv (+ zonal + ensemble) -> CF-NetCDF for presto-viz.

The Temperature 12k reconstruction is ZONAL (six 30-deg latitude bands). To
drive presto-viz (which expects a spatial field) we emit a latitudinal-stripe
pseudo-grid: each band's consensus median fills its latitudes across all
longitudes. visualize.yml then autodetects lat+lon dims and routes to presto-viz.

Variables (matching presto-viz's tas / tas_gm contract):
  tas    (time, lat, lon)  -- zonal-stripe consensus median field (degC)
  tas_gm (time, ens)       -- pooled multi-method global-mean ensemble (degC)

If the zonal/ensemble files are absent (e.g. a single-method run), falls back to
a 1D (time-only) NetCDF -> static-Pages tile UI.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

N_BANDS = 6
LAT = np.arange(-87.5, 90.0, 5.0)     # 36 cells, 5-deg centers
LON = np.arange(2.5, 360.0, 5.0)      # 72 cells, 5-deg centers


def _band_of(lat):
    return min(int((lat + 90) // 30) + 1, N_BANDS)   # 1..6


def convert(in_csv: Path, out_nc: Path) -> None:
    d = pd.read_csv(in_csv)
    zonal_p = in_csv.parent / "reconstruction_zonal.csv"
    ens_p = in_csv.parent / "reconstruction_ensemble.csv"

    if zonal_p.exists() and ens_p.exists():
        _write_spatial(d, zonal_p, ens_p, out_nc)
    else:
        _write_1d(d, out_nc)


def _write_spatial(d, zonal_p, ens_p, out_nc):
    z = pd.read_csv(zonal_p)
    e = pd.read_csv(ens_p)
    age = z["binAges"].to_numpy(float)
    year = (1950 - age).astype(int)
    order = np.argsort(year)
    year = year[order]

    # zonal-stripe spatial field: tas[t, lat, lon] = band_median[t, band(lat)]
    band_med = {b: z[f"band{b}_med"].to_numpy(float)[order] for b in range(1, N_BANDS + 1)}
    tas = np.full((year.size, LAT.size, LON.size), np.nan)
    for li, lat in enumerate(LAT):
        tas[:, li, :] = band_med[_band_of(lat)][:, None]

    # pooled global-mean ensemble -> (time, ens)
    gm = e.drop(columns=["binAges"]).to_numpy(float)[order]   # (time, ens)

    ds = xr.Dataset(
        data_vars={
            "tas":    (("time", "lat", "lon"), tas),
            "tas_gm": (("time", "ens"), gm),
        },
        coords={"time": year, "lat": LAT, "lon": LON, "ens": np.arange(gm.shape[1]) + 1},
        attrs={
            "title": "Temperature 12k multi-method composite (zonal stripes)",
            "source": "presto-Temp12k_Composites (SCC/DCC/GAM/CPS/PaiCo consensus)",
            "reference": "Kaufman et al. 2020, Sci. Data 7:201, doi:10.1038/s41597-020-0530-7",
            "note": "Spatial field is a 30-deg zonal-band reconstruction shown as latitudinal stripes.",
            "Conventions": "CF-1.10",
        },
    )
    ds["time"].attrs.update(units="years_AD", standard_name="time", long_name="year CE")
    ds["lat"].attrs.update(units="degrees_north", standard_name="latitude")
    ds["lon"].attrs.update(units="degrees_east", standard_name="longitude")
    ds["tas"].attrs.update(long_name="zonal-band temperature anomaly", units="degC")
    ds["tas_gm"].attrs.update(long_name="global-mean temperature anomaly (ensemble)", units="degC")
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_nc)
    print(f"[outputs_to_netcdf] wrote SPATIAL {out_nc} "
          f"(time={year.size}, lat={LAT.size}, lon={LON.size}, ens={gm.shape[1]})")


def _write_1d(d, out_nc):
    d = d.sort_values("year").reset_index(drop=True)

    def col(name):
        return d[name].to_numpy(float) if name in d.columns else np.full(len(d), np.nan)

    ds = xr.Dataset(
        data_vars={"gmst": ("time", col("mean")),
                   "gmst_lo_95": ("time", col("lo_95")),
                   "gmst_hi_95": ("time", col("hi_95"))},
        coords={"time": d["year"].to_numpy(int)},
        attrs={"title": "Temperature 12k composite GMST (1D)", "Conventions": "CF-1.10"},
    )
    ds["time"].attrs.update(units="years_AD", standard_name="time")
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_nc)
    print(f"[outputs_to_netcdf] wrote 1D {out_nc} ({len(d)} time steps)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", required=True, type=Path)
    ap.add_argument("--out-nc", required=True, type=Path)
    args = ap.parse_args()
    convert(args.in_csv, args.out_nc)


if __name__ == "__main__":
    main()
