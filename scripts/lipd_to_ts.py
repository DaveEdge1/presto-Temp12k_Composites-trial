#!/usr/bin/env python3
"""Temp12k LiPD legacy pickle -> proxy_ts.json (full-fidelity timeseries).

The five Temperature-12k composite methods (SCC, DCC, GAM, CPS, PaiCo) all
operate on individual proxy time series with their native irregular age axes
and per-record interpretation metadata -- NOT on a year x proxy matrix. So,
unlike the presto-template demo adapter (which flattens to an annual matrix
and loses ages + interpretation), this adapter preserves every record's full
`age`/`value` arrays and the metadata the composite engines need.

Input: the lipdverse "legacy" pickle, e.g. Temp12k1_0_2.pkl. It carries both
a `D` dict of datasets and a pre-flattened `TS` list (one entry per paleo
column). We use `TS` directly -- it is exactly what the original R analysis
consumed via `lipdR::extractTs` + `pullTsVariable`.

Output: proxy_ts.json -- a list of record objects, each:
    {
      "id": <unique>, "dataSetName", "variableName",
      "age": [...yr BP...], "values": [...],
      "lat", "lon", "elev", "archiveType",
      "units", "inCompilation",
      "proxy", "proxyGeneral", "proxyType",
      "direction", "scope", "seasonality", "seasonalityGeneral",
      "variable", "variableGroup",
      "uncertainty1sd": <float or null>
    }

Records kept: those in the Temp12k compilation with a temperature
interpretation. Method-specific subsetting (degC-only for SCC/DCC; seasonality
{annual,summerOnly,winterOnly}) is applied downstream so the JSON stays a
faithful superset.
"""

from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
from pathlib import Path

try:
    import yaml  # only needed when --uncertainties is supplied
except Exception:  # pragma: no cover
    yaml = None


def _num_list(x):
    """Coerce a sequence to a list of floats (NaN for non-numeric)."""
    if x is None:
        return None
    out = []
    try:
        for v in x:
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                out.append(float("nan"))
    except TypeError:
        return None
    return out


def _safe_float(x):
    try:
        f = float(x)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _interp0(rec: dict) -> dict:
    it = rec.get("paleoData_interpretation")
    if isinstance(it, list) and it and isinstance(it[0], dict):
        return it[0]
    if isinstance(it, dict):
        return it
    return {}


def _coerce_age(rec: dict):
    """Return age as a list in yr BP (1950 reference).

    Temp12k records carry `age` in yr BP (ageUnits 'BP'). A few may instead
    carry `year` (yr AD/CE) -> convert to BP.
    """
    age = _num_list(rec.get("age"))
    if age is not None and any(math.isfinite(v) for v in age):
        units = str(rec.get("ageUnits", "") or "").lower()
        if "ad" in units or "ce" in units:  # age stored as calendar year
            return [1950.0 - v for v in age]
        return age
    year = _num_list(rec.get("year"))
    if year is not None and any(math.isfinite(v) for v in year):
        return [1950.0 - v for v in year]
    return None


def _is_temperature(rec: dict, interp: dict) -> bool:
    if str(interp.get("variableGroup", "")).lower() == "temperature":
        return True
    if str(interp.get("variable", "")).strip().upper() == "T":
        return True
    if str(rec.get("paleoData_units", "")).lower() == "degc":
        return True
    return False


def _in_temp12k(rec: dict) -> bool:
    a = str(rec.get("paleoData_inCompilation", "") or "").lower()
    b = str(rec.get("paleoData_inCompilationBeta", "") or "").lower()
    return "temp12k" in a or "temp12k" in b


def _proxy_type(rec: dict) -> str:
    archive = str(rec.get("archiveType", "unknown") or "unknown").lower()
    proxy = str(rec.get("paleoData_proxyGeneral")
                or rec.get("paleoData_proxy") or "").lower()
    return f"{archive}.{proxy}" if proxy else archive


def _record_uncertainty(rec: dict):
    """Per-record 1-sigma temperature uncertainty if the record states one."""
    for k in ("paleoData_temperature12kUncertainty",
              "paleoData_calibration_uncertainty",
              "paleoData_calibration_calibrationUncertainty"):
        v = _safe_float(rec.get(k))
        if v is not None and v > 0:
            return v
    return None


def _load_unc_table(path: Path):
    if path is None:
        return None
    if yaml is None:
        print("[lipd_to_ts] pyyaml missing; cannot read uncertainties table",
              file=sys.stderr)
        return None
    data = yaml.safe_load(path.read_text()) or {}
    # normalize keys to lowercase for matching
    return {str(k).lower(): v for k, v in data.items()}


# Map the pickle's proxy / proxyGeneral spellings onto the Table-2 (Summary
# sheet) proxy keys in proxy_uncertainties.yml.
_PROXY_ALIASES = {
    "mg/ca": "mgca",
    "mgca": "mgca",
    "dinocyst": "other microfossils/dinocyst",
    "foraminifera": "other microfossils/foraminifera",
    "diatom": "other microfossils/diatoms",
    "diatoms": "other microfossils/diatoms",
    "radiolaria": "other microfossils/radiolaria",
    "isotope": "d18o",
    "d18o": "d18o",
    "tex86": "gdgt (tex86)",
    "gdgt": "gdgt (mbt/cbt as well as brgdgt fractional abundance)",
    "uk37": "alkenone",
    "alkenone": "alkenone",
    "pollen": "pollen",
    "chironomid": "chironomid",
}


def _season_bucket(season_general: str) -> str:
    s = (season_general or "").lower()
    if s.startswith("summer"):
        return "summer"
    if s.startswith("winter"):
        return "winter"
    return "annual"


def _lookup_unc(unc_table, archive: str, proxy: str, season: str, default: float):
    """Table-2 style lookup: try archive.proxy, then proxy (with aliases), else default."""
    if not unc_table:
        return default
    archive = (archive or "").lower()
    proxy = (proxy or "").lower()
    season = _season_bucket(season)
    candidates = [f"{archive}.{proxy}", proxy, _PROXY_ALIASES.get(proxy, proxy), archive]
    for key in candidates:
        entry = unc_table.get(key)
        if entry is None:
            continue
        if isinstance(entry, dict):
            for sk in (str(season).lower(), "annual", "value"):
                if sk in entry and _safe_float(entry[sk]) is not None:
                    return float(entry[sk])
        else:
            v = _safe_float(entry)
            if v is not None:
                return v
    return default


def build(pkl_path: Path, out_json: Path, unc_path: Path = None,
          default_unc: float = 2.1) -> None:
    with pkl_path.open("rb") as f:
        raw = pickle.load(f)

    TS = raw.get("TS") if isinstance(raw, dict) else None
    if not TS:
        raise SystemExit("Pickle has no 'TS' list; expected a lipdverse legacy pickle.")

    unc_table = _load_unc_table(unc_path) if unc_path else None

    out = []
    seen = set()
    n_total = len(TS)
    n_not_temp12k = n_not_temp = n_bad_series = 0

    for rec in TS:
        if not _in_temp12k(rec):
            n_not_temp12k += 1
            continue
        interp = _interp0(rec)
        if not _is_temperature(rec, interp):
            n_not_temp += 1
            continue

        values = _num_list(rec.get("paleoData_values"))
        age = _coerce_age(rec)
        if values is None or age is None or len(values) != len(age):
            n_bad_series += 1
            continue
        if not any(math.isfinite(v) for v in values):
            n_bad_series += 1
            continue

        rid = str(rec.get("paleoData_TSid")
                  or f'{rec.get("dataSetName")}__{rec.get("paleoData_variableName")}')
        base, i = rid, 1
        while rid in seen:
            i += 1
            rid = f"{base}__{i}"
        seen.add(rid)

        archive = str(rec.get("archiveType", "unknown") or "unknown")
        proxy = str(rec.get("paleoData_proxyGeneral")
                    or rec.get("paleoData_proxy") or "")
        season_general = str(interp.get("seasonalityGeneral", "") or "")

        unc = _record_uncertainty(rec)
        if unc is None:
            unc = _lookup_unc(unc_table, archive, proxy, season_general, default_unc)

        # JSON has no NaN literal that R's jsonlite will read -> emit null,
        # which jsonlite parses back to NA in a numeric vector.
        age_j = [None if (v is None or not math.isfinite(v)) else v for v in age]
        values_j = [None if (v is None or not math.isfinite(v)) else v for v in values]

        out.append({
            "id": rid,
            "dataSetName": str(rec.get("dataSetName", "") or ""),
            "variableName": str(rec.get("paleoData_variableName", "") or ""),
            "age": age_j,
            "values": values_j,
            "lat": _safe_float(rec.get("geo_meanLat")),
            "lon": _safe_float(rec.get("geo_meanLon")),
            "elev": _safe_float(rec.get("geo_meanElev")),
            "archiveType": archive,
            "units": str(rec.get("paleoData_units", "") or ""),
            "inCompilation": str(rec.get("paleoData_inCompilation")
                                 or rec.get("paleoData_inCompilationBeta") or ""),
            "proxy": str(rec.get("paleoData_proxy", "") or ""),
            "proxyGeneral": str(rec.get("paleoData_proxyGeneral", "") or ""),
            "proxyType": _proxy_type(rec),
            "direction": str(interp.get("direction", "") or ""),
            "scope": str(interp.get("scope", "") or ""),
            "seasonality": str(interp.get("seasonality", "") or ""),
            "seasonalityGeneral": season_general,
            "variable": str(interp.get("variable", "") or ""),
            "variableGroup": str(interp.get("variableGroup", "") or ""),
            "uncertainty1sd": unc,
        })

    if not out:
        raise SystemExit("No Temp12k temperature records extracted -- cannot proceed.")

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, allow_nan=False))

    # diagnostics (stderr)
    n_degc = sum(1 for r in out if r["units"].lower() == "degc")
    seasons = {}
    for r in out:
        seasons[r["seasonalityGeneral"]] = seasons.get(r["seasonalityGeneral"], 0) + 1
    print(f"[lipd_to_ts] {n_total} TS rows -> {len(out)} Temp12k temperature records "
          f"({n_degc} degC). dropped: not-temp12k={n_not_temp12k}, "
          f"not-temperature={n_not_temp}, bad-series={n_bad_series}", file=sys.stderr)
    print(f"[lipd_to_ts] seasonalityGeneral counts: {seasons}", file=sys.stderr)
    print(f"[lipd_to_ts] wrote {out_json}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickle", required=True, type=Path)
    ap.add_argument("--out-json", required=True, type=Path)
    ap.add_argument("--uncertainties", type=Path, default=None,
                    help="reference_data/proxy_uncertainties.yml (Table 2)")
    ap.add_argument("--default-unc", type=float, default=2.1,
                    help="fallback 1-sigma uncertainty (degC) when no match")
    args = ap.parse_args()
    build(args.pickle, args.out_json, args.uncertainties, args.default_unc)


if __name__ == "__main__":
    main()
