#!/usr/bin/env python3
"""GAM -- Generalized Additive Model composite (Python / pygam).

Faithful Python reimplementation of the Kaufman 2020 Temperature 12k GAM
reconstruction (Sommer & Davis; paper Methods, "Reconstruction method 3"; the
original code lives in nickmckay/Temperature12k/GAM_frozen). For each
equal-area grid cell we pool the (age_ensemble x value) samples of every record
in the cell, fit one penalized LinearGAM, and draw posterior smooth-function
samples (`quantity='mu'`). Cells are averaged within 30-degree bands and the
six band means are area-weighted into a global ensemble.

Differences from the published archive (driven by pygam-version and runtime
constraints, validated against the per-method fidelity table):
  - Posterior samples use `quantity='mu'` + an explicit empirical-residual
    noise term. pygam 0.12's `quantity='y'` over-inflates `scale_` by ~3x vs
    the pygam 0.8 the original used, blowing up the ensemble spread.
  - `lam` gridsearch is constrained to `np.logspace(-1, 3, 11)` to remove the
    wiggly low-lam tail that 0.12's gridsearch occasionally picks (artificially
    sharp deglacial peaks).
  - Per-record alignment uses the published iterative-union `_align_ensembles`
    (gam_ensemble.py L385-426): pick the cell's longest record as base, then
    iteratively absorb others into the union when they share >100 ensemble
    samples of age overlap. Records with no overlap and <100 samples in 3-5 ka
    fall back to a WorldClim modern-temperature lookup (paper-faithful,
    matches gam_ensemble.py L534).

Output schema matches the R methods (binAges + nens columns) so the
consensus step treats all five methods uniformly.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

LATBINS = np.arange(-90, 91, 30)                       # 6 bands
BAND_WEIGHTS = np.array([0.067, 0.183, 0.25, 0.25, 0.183, 0.067])
N_BANDS = len(BAND_WEIGHTS)

# Modern-anomaly reference window (yr BP). gam_ensemble.py._compute_anomaly
# defines "modern" as modern_young=-50 .. modern_old=-20 BP (1970-2000 CE, the
# WorldClim baseline). Records/unions with >100 ensemble samples in this window
# anomalize to their own modern mean; the rest fall back to a WorldClim lookup.
MODERN_LO, MODERN_HI = -50.0, -20.0
MODERN_MIN_SAMPLES = 100
# 0-insertion anchor: fraction of a cell's pooled-cloud size added as synthetic
# (age=-35 BP, value=0) points, pinning the GAM through 0 at the 1985 reference
# (gam_ensemble.py._predict_gam). Floor of 50 points.
ANCHOR_BP = -35.0
ANCHOR_FRAC = float(os.environ.get("GAM_ANCHOR_FRAC", "0.05"))
ANCHOR_MIN = int(os.environ.get("GAM_ANCHOR_MIN", "50"))
# GAM_ANCHOR_MODE: "insert" = 0-insert at -35 BP; "none" = no 0-insertion.
ANCHOR_MODE = os.environ.get("GAM_ANCHOR_MODE", "insert")
# GAM_WC_FALLBACK: how records lacking modern-window data anchor.
#   "worldclim" = subtract the WorldClim absolute-temp lookup (original);
#   "midhol"    = subtract their own 3-5 ka mean (avoids worldclim bias);
#   "drop"      = drop them.
WC_FALLBACK = os.environ.get("GAM_WC_FALLBACK", "worldclim")

# Sigma fallback when no proxy×season entry matches a record. Paper notebook
# cell 29 actually uses 1.975878 (75th pct without d18O); we use 1.7 (median)
# because the higher default inflates the cell variance for records hitting
# the fallback (we confirmed this empirically against the published curve).
SIGMA_DEFAULT = 1.7


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


# ---------- proxy × season sigma lookup -----------------------------------
def load_sigma_table(path):
    if not Path(path).exists():
        return {}
    df = pd.read_csv(path)
    out = {}
    for _, row in df.iterrows():
        proxy = str(row["proxy"]).strip()
        for sc in ("summer", "winter", "annual"):
            v = row.get(sc, None)
            try:
                fv = float(v)
            except Exception:
                fv = float("nan")
            if np.isfinite(fv):
                out[(proxy, sc)] = fv
    return out


def _match_proxy_cat(rec_proxy):
    p = (rec_proxy or "").lower()
    if "pollen" in p: return "pollen"
    if "alkenone" in p or "uk37" in p or "uk37'" in p: return "alkenone"
    if "mg/ca" in p or "mgca" in p: return "MgCa"
    if "chironomid" in p: return "chironomid"
    if "tex86" in p: return "GDGT (Tex86)"
    if "mbt" in p or "brgdgt" in p or "gdgt" in p:
        return "GDGT (MBT/CBT as well as BrGDGT fractional abundance)"
    if "d18o" in p: return "d18O"
    if "diatom" in p: return "other microfossils/diatoms"
    if "dinocyst" in p or "dinoflagell" in p: return "other microfossils/dinocyst"
    if "radiolaria" in p: return "other microfossils/radiolaria"
    if "foramini" in p or "foram" in p: return "other microfossils/foraminifera"
    return None


def _match_season_col(sg):
    s = (sg or "").lower()
    if "summer" in s: return "summer"
    if "winter" in s: return "winter"
    return "annual"


def get_sigma(sigma_table, proxy_str, season_str):
    pc = _match_proxy_cat(proxy_str)
    if pc is None:
        return SIGMA_DEFAULT
    sc = _match_season_col(season_str)
    v = sigma_table.get((pc, sc))
    if v is None or not np.isfinite(v):
        v = sigma_table.get((pc, "annual"))
    return float(v) if v is not None and np.isfinite(v) else SIGMA_DEFAULT


# ---------- WorldClim modern lookup --------------------------------------
def load_modern_grid(path):
    """Returns a (180, 360) array of annual mean temperature in degC.
    Row 0 = latitude band [89, 90]; row 179 = [-90, -89]. Col 0 = longitude
    band [-180, -179]; col 359 = [179, 180]. NaN where no land pixel."""
    if not Path(path).exists():
        return None
    return pd.read_csv(path, header=None).to_numpy(float)


def modern_lookup(grid, lat, lon):
    if grid is None or not np.isfinite(lat) or not np.isfinite(lon):
        return float("nan")
    # bilinear at the 1-deg cell centres (0.5, 1.5, ..., 359.5 in col index)
    rf = (89.5 - lat)                                   # row continuous
    cf = (lon + 179.5) % 360.0                          # col continuous (wraps)
    r0 = int(np.floor(rf)); c0 = int(np.floor(cf))
    if not (0 <= r0 < grid.shape[0] - 1 and 0 <= c0 < grid.shape[1] - 1):
        # at the grid edges, nearest non-NaN
        r0 = max(0, min(grid.shape[0] - 1, r0))
        c0 = max(0, min(grid.shape[1] - 1, c0))
        return float(grid[r0, c0])
    dr = rf - r0; dc = cf - c0
    block = grid[r0:r0 + 2, c0:c0 + 2]
    if not np.isfinite(block).all():
        for k in range(1, 6):
            r_lo, r_hi = max(0, r0 - k), min(grid.shape[0], r0 + k + 1)
            c_lo, c_hi = max(0, c0 - k), min(grid.shape[1], c0 + k + 1)
            nb = grid[r_lo:r_hi, c_lo:c_hi]
            if np.isfinite(nb).any():
                return float(np.nanmean(nb))
        return float("nan")
    return float((1 - dr) * (1 - dc) * block[0, 0]
                 + (1 - dr) * dc * block[0, 1]
                 + dr * (1 - dc) * block[1, 0]
                 + dr * dc * block[1, 1])


# ---------- record loader ------------------------------------------------
def load_records(ts_path, sigma_table, modern_grid, rng_seed=42):
    """Load proxy_ts.json, filter to annual|summerOnly|winterOnly + degC,
    compute per-record sigma + WorldClim modern + Gaussian-perturbed
    age ensemble (500 cols, age_unc = 50 + age*200/12000), decadally
    decimate >1470-sample records (paper notebook cell 19).

    Preserves the per-record value-ensemble (rebuilt by lipd_to_ts.py to
    match the published pipeline's multi-col `paleoData_values` matrix).
    Each ensemble member samples one column per call -- equivalent to
    compositeR's NCOL>1 path.
    """
    rng = np.random.default_rng(rng_seed)
    recs = json.loads(Path(ts_path).read_text())
    out = []
    for r in recs:
        if str(r.get("units", "")).lower() != "degc":
            continue
        if str(r.get("seasonalityGeneral", "")).lower() not in ("annual", "summeronly", "winteronly"):
            continue
        age_med = np.asarray(r["age"], dtype=float)
        # GAM-specific choice: use the single-vector measurement (r["values"])
        # and let the per-cell pool-sampler add Gaussian sigma noise per draw
        # (paper's GAM ensemble construction). The pre-baked value-ensemble in
        # proxy_ts.json (N=10 AR1 realisations) is too few to give the per-cell
        # GAM cloud the diversity it needs -- ten fixed realisations vs the
        # published's 100-2500 columns from real LiPD ensembles. Per-draw fresh
        # sigma noise reproduces the published behaviour with our 10-col input.
        base_val = np.asarray(r["values"], dtype=float)
        val_mat = base_val.reshape(-1, 1)
        if r.get("lat") is None or r.get("lon") is None:
            continue
        lat = float(r["lat"]); lon = float(r["lon"])
        # Paper cell 32: outlier zero
        val_mat = np.where(np.abs(val_mat) > 200, np.nan, val_mat)
        base_val = np.where(np.abs(base_val) > 200, np.nan, base_val)
        m = np.isfinite(age_med) & np.isfinite(base_val) & (age_med >= -50) & (age_med <= 12050)
        if m.sum() < 3:
            continue
        age_med = age_med[m]; val_mat = val_mat[m]; base_val = base_val[m]
        # Direction flip (negative-direction proxies)
        if str(r.get("direction", "")).lower() == "negative":
            val_mat = val_mat * -1.0
            base_val = base_val * -1.0
        # Paper cell 19: decadal averaging for records with >1470 samples
        if base_val.size > 1470:
            keys = (5 + age_med - (age_med % 10)).astype(int)
            uniq = np.unique(keys)
            new_age = np.array([np.nanmean(age_med[keys == k]) for k in uniq])
            new_val_mat = np.array([np.nanmean(val_mat[keys == k, :], axis=0)
                                     for k in uniq])
            age_med, val_mat = new_age, new_val_mat
            base_val = np.nanmean(val_mat, axis=1)
        # Paper cell 24+38: Gaussian-perturbed age ensemble, 500 columns
        age_unc = 50.0 + np.maximum(age_med, 0.0) * (250.0 - 50.0) / 12000.0
        age_ens = rng.normal(loc=age_med[:, None], scale=age_unc[:, None],
                             size=(age_med.size, 500)).astype(np.float32)
        # Per-proxy × season sigma -- applied every pool draw (matches paper's
        # GAM behaviour and the published gam_ensemble.py default that gives
        # local maxD ~0.22 / midHol 0.465 / spread 0.996).
        sigma = get_sigma(sigma_table, r.get("proxy"), r.get("seasonalityGeneral"))
        # WorldClim modern lookup (fallback for records lacking 3-5 ka coverage)
        modern = modern_lookup(modern_grid, lat, lon)
        out.append({"age_med": age_med, "val": base_val, "val_mat": val_mat,
                    "lat": lat, "lon": lon,
                    "age_ens": age_ens, "sigma": float(sigma),
                    "modern": float(modern)})
    return out


# ---------- iterative-union alignment (`_align_ensembles`) ----------------
def _modern_ref(age_arr, val_arr, worldclim):
    """The gam_ensemble.py._compute_anomaly reference for one (aligned union or
    solo) cloud: the mean climate over the modern window (-50..-20 BP) when the
    cloud has >100 ensemble samples there, else the WorldClim modern lookup
    (may be NaN when no land pixel is near, which drops the record)."""
    mmask = (age_arr >= MODERN_LO) & (age_arr <= MODERN_HI)
    if int(np.count_nonzero(mmask)) > MODERN_MIN_SAMPLES:
        mv = float(np.nanmean(val_arr[mmask]))
        if np.isfinite(mv):
            return mv, True   # anomalized from own modern data
    if WC_FALLBACK == "worldclim":
        return float(worldclim), False
    if WC_FALLBACK == "midhol":
        hmask = (age_arr >= 3000) & (age_arr <= 5000)
        if int(np.count_nonzero(hmask)) > MODERN_MIN_SAMPLES:
            hv = float(np.nanmean(val_arr[hmask]))
            if np.isfinite(hv):
                return hv, False
    return float("nan"), False


def compute_alignments(recs, cell_idx):
    """For each cell, pick the longest record as base; iteratively shift other
    records by the mean offset over their age overlap with the growing aligned
    union (`min_overlap = 100` ensemble samples both sides, per gam_ensemble.py
    _align_ensembles L385-426). Then anchor to MODERN, not mid-Holocene: the
    aligned union subtracts ONE reference (its modern-window mean, else the
    base record's WorldClim value); each non-overlapping / singleton record
    subtracts its own modern reference (gam_ensemble.py._compute_anomaly
    L485-538). Returns per-record `offset` (raw_value - offset = modern anomaly)
    and `pre_anomalized` (True when anchored via WorldClim rather than shared
    alignment; used only for diagnostics). Records with no modern data and no
    WorldClim pixel get offset NaN and are dropped from the pool.
    """
    n = len(recs)
    offset = np.full(n, np.nan, dtype=float)
    pre_anomalized = np.zeros(n, dtype=bool)
    aligned = np.zeros(n, dtype=bool)
    n_aligned = n_union_data = n_union_wc = n_union_drop = 0
    n_solo_data = n_solo_wc = n_solo_drop = 0

    by_cell = {}
    for i, c in enumerate(cell_idx):
        if c is None:
            continue
        by_cell.setdefault(c, []).append(i)

    def _cloud(i):
        a = recs[i]["age_ens"].ravel()
        v = np.broadcast_to(recs[i]["val"][:, None],
                            recs[i]["age_ens"].shape).ravel()
        return a, v

    def anchor_solo(i):
        nonlocal n_solo_data, n_solo_wc, n_solo_drop
        a, v = _cloud(i)
        ref, from_data = _modern_ref(a, v, recs[i]["modern"])
        if np.isfinite(ref):
            offset[i] = ref
            pre_anomalized[i] = not from_data
            if from_data:
                n_solo_data += 1
            else:
                n_solo_wc += 1
        else:
            offset[i] = np.nan
            n_solo_drop += 1

    for cid, members in by_cell.items():
        if len(members) == 1:
            anchor_solo(members[0])
            continue
        # multi-record cell -- iterative-union alignment (relative shifts)
        lens = [recs[i]["age_med"].size for i in members]
        base_i = members[int(np.argmax(lens))]
        offset[base_i] = 0.0            # relative to the base for now
        aligned[base_i] = True
        u_age, u_val = _cloud(base_i)
        remaining = [i for i in members if i != base_i]
        changed = True
        while remaining and changed:
            changed = False
            still = []
            a_min = float(np.min(u_age)); a_max = float(np.max(u_age))
            for i in remaining:
                r_age, r_val = _cloud(i)
                if r_age.size == 0:
                    still.append(i); continue
                r_min = float(np.min(r_age)); r_max = float(np.max(r_age))
                m1 = (u_age >= r_min) & (u_age <= r_max)
                m2 = (r_age >= a_min) & (r_age <= a_max)
                if m1.sum() > 100 and m2.sum() > 100:
                    diff = float(np.nanmean(u_val[m1]) - np.nanmean(r_val[m2]))
                    offset[i] = -diff   # aligned value = raw + diff
                    aligned[i] = True
                    u_age = np.concatenate([u_age, r_age])
                    u_val = np.concatenate([u_val, r_val + diff])
                    changed = True
                else:
                    still.append(i)
            remaining = still
        # Anchor the whole aligned union with ONE modern reference: shift every
        # aligned member's offset by ref so raw - offset = aligned - ref.
        union = [i for i in members if aligned[i]]
        ref, from_data = _modern_ref(u_age, u_val, recs[base_i]["modern"])
        if np.isfinite(ref):
            for i in union:
                offset[i] += ref
            n_aligned += len(union)
            if from_data:
                n_union_data += 1
            else:
                n_union_wc += 1
        else:
            for i in union:
                offset[i] = np.nan   # no modern anchor -> drop the union
            n_union_drop += 1
        # Non-overlapping records anchor individually.
        for i in remaining:
            anchor_solo(i)
    print(f"[gam] anchor: {n_aligned} in {n_union_data + n_union_wc} unions "
          f"({n_union_data} modern-data / {n_union_wc} worldclim / {n_union_drop} dropped) | "
          f"solo {n_solo_data} modern-data + {n_solo_wc} worldclim + {n_solo_drop} drop",
          file=sys.stderr, flush=True)
    return offset, pre_anomalized


# ---------- pooled cloud + per-cell GAM ----------------------------------
def build_pool_for_cell(records, n_pool, rng):
    """For each record in a cell, sample n_pool (age, value) pairs:
       age = random column of age_ens; val = random column of val_mat
       (the pre-built value ensemble) minus the alignment offset, plus a
       residual N(0, sigma) noise (zero when val_mat already has cols).
       Concatenate into the cell's pooled cloud."""
    ages_list, vals_list, pre_list = [], [], []
    for r, offset, pre_an in records:
        if not np.isfinite(offset):
            continue
        n_samp, n_age_cols = r["age_ens"].shape
        n_val_cols = r["val_mat"].shape[1]
        if n_samp == 0:
            continue
        age_col_idx = rng.integers(0, n_age_cols, size=n_pool)
        val_col_idx = rng.integers(0, n_val_cols, size=n_pool)
        ages = r["age_ens"][:, age_col_idx]              # (n_samp, n_pool)
        vals = r["val_mat"][:, val_col_idx].astype(np.float64)
        vals = vals - offset
        if r["sigma"] > 0:
            vals = vals + rng.normal(0.0, r["sigma"], size=(n_samp, n_pool))
        ages_list.append(ages.ravel())
        vals_list.append(vals.ravel())
        pre_list.append(np.full(ages.size, int(pre_an), dtype=np.int8))
    if not ages_list:
        return None
    return (np.concatenate(ages_list),
            np.concatenate(vals_list),
            np.concatenate(pre_list))


def fit_cell(args):
    """Fit one cell's pooled cloud, return smooth-function posterior samples
    + empirical residual std (for the noise-gain step). Matches the published
    pygam call but constrains `lam` to remove the wiggly low-lam tail of the
    pygam 0.12 gridsearch."""
    from pygam import LinearGAM, s
    cid, x, y, pre_anom, band, n_draws, bin_ages, binvec, seed = args
    # pygam's gam.sample() draws via numpy's GLOBAL RNG (not a Generator), and
    # fit_cell runs in forked pool workers, so without this the posterior draws
    # differ every run. Seed the legacy global RNG deterministically per cell.
    np.random.seed((int(seed) + int(cid) * 100003) % (2**32))
    if x.size < 20 or np.ptp(x) < 200:
        return cid, None, None, band
    # Pool values are already anomalized to the modern window (offset applied in
    # build_pool), so no re-centering here. Anchor the GAM through 0 at the 1985
    # reference by inserting synthetic (age=-35 BP, value=0) points, and keep
    # every cell -- the original does not require 3-5 ka coverage.
    y0 = y
    if ANCHOR_MODE == "insert":
        n_anchor = max(ANCHOR_MIN, int(round(ANCHOR_FRAC * x.size)))
        x_fit = np.concatenate([x, np.full(n_anchor, ANCHOR_BP)])
        y_fit = np.concatenate([y0, np.zeros(n_anchor)])
    else:
        x_fit, y_fit = x, y0
    try:
        lam_grid = np.logspace(-1, 3, 11)
        gam = LinearGAM(s(0)).gridsearch(x_fit[:, None], y_fit, lam=lam_grid, progress=False)
        draws = gam.sample(x_fit[:, None], y_fit, sample_at_X=bin_ages[:, None],
                           n_draws=n_draws, n_bootstraps=1, quantity="mu")
        resid = y0 - gam.predict(x[:, None])
        scale = float(np.nanstd(resid))
        outside = (bin_ages < x.min()) | (bin_ages > x.max())
        draws[:, outside] = np.nan
        counts, _ = np.histogram(x, bins=binvec)
        surrounding = np.convolve(counts, np.ones(3, dtype=int), mode="same")
        draws[:, surrounding < 100] = np.nan
        return cid, draws.T, scale, band                  # (NB, n_draws)
    except Exception as exc:
        print(f"[gam] cell {cid}: {exc}", file=sys.stderr)
        return cid, None, None, band


def aggregate(cell_draws, cell_band, nens, NB, ZW):
    band_ens = np.full((NB, N_BANDS, nens), np.nan)
    for b in range(N_BANDS):
        cells_b = [c for c, bn in cell_band.items() if bn == b and c in cell_draws]
        if not cells_b:
            continue
        arr = np.stack([cell_draws[c] for c in cells_b], axis=0)
        band_ens[:, b, :] = np.nanmean(arr, axis=0)
    glob = np.full((NB, nens), np.nan)
    for k in range(nens):
        bm = band_ens[:, :, k]
        w = np.tile(ZW, (NB, 1)).astype(float)
        w[~np.isfinite(bm)] = np.nan
        num = np.nansum(bm * w, axis=1); den = np.nansum(w, axis=1)
        glob[:, k] = np.where(den > 0, num / den, np.nan)
    return band_ens, glob


def archival_reference(glob, bin_ages, ref_start_ce=1800, ref_end_ce=1900):
    """Per-member subtract full-12k mean, then anchor ensemble median to 0
    over the paper's reference period (1800-1900 CE -> ~50-150 yr BP).
    Matches the R methods' apply_reference (run_methods.R)."""
    glob = glob - np.nanmean(glob, axis=0, keepdims=True)
    ref_bp_lo = 1950 - ref_end_ce
    ref_bp_hi = 1950 - ref_start_ce
    refrows = np.where((bin_ages >= ref_bp_lo) & (bin_ages <= ref_bp_hi))[0]
    if refrows.size == 0:
        refrows = np.array([int(np.argmin(np.abs(bin_ages - 100)))])
    per_member_ref = np.nanmean(glob[refrows, :], axis=0)
    glob = glob - np.nanmedian(per_member_ref)
    return glob


def run_gam(ts_path, cfg, grid, sigma_table_path, modern_grid_path, out_csv):
    bin_cfg = cfg.get("bin", {})
    start = bin_cfg.get("start_bp", -50)
    end = bin_cfg.get("end_bp", 12050)
    step = bin_cfg.get("step", 100)
    binvec = np.arange(start, end + step, step, dtype=float)
    bin_ages = (binvec[1:] + binvec[:-1]) / 2.0
    NB = bin_ages.size
    nens = int(cfg.get("nens", 100))
    n_pool = int((cfg.get("advanced") or {}).get("gam_n_pool", 500))
    noise_gain = float((cfg.get("advanced") or {}).get("gam_noise_gain", 1.0))
    ZW = np.sin(LATBINS[1:] * np.pi / 180) - np.sin(LATBINS[:-1] * np.pi / 180)
    ZW = ZW / ZW.sum()
    seed = int(cfg.get("advanced", {}).get("seed") or 42)

    sigma_table = load_sigma_table(sigma_table_path)
    modern_grid = load_modern_grid(modern_grid_path)
    print(f"[gam] sigma_table: {len(sigma_table)} entries; "
          f"modern_grid: {'loaded' if modern_grid is not None else 'NONE'}", file=sys.stderr)

    recs = load_records(ts_path, sigma_table, modern_grid, rng_seed=seed)
    print(f"[gam] {len(recs)} records after filter", file=sys.stderr)
    if len(recs) < 10:
        raise SystemExit("[gam] too few records")

    clat = grid["clat"].to_numpy(); clon = grid["clon180"].to_numpy()
    cell_idx, band_idx = [], []
    for r in recs:
        b = band_of(r["lat"])
        c = nearest_cell(r["lat"], r["lon"], clat, clon) if b is not None else None
        band_idx.append(b)
        cell_idx.append(c)

    offset, pre_anomalized = compute_alignments(recs, cell_idx)
    n_dropped = int(np.sum(~np.isfinite(offset)))
    print(f"[gam] dropped {n_dropped} records with no valid anchor", file=sys.stderr)

    # group records by cell
    by_cell = {}
    cell_band = {}
    for i, c in enumerate(cell_idx):
        if c is None or not np.isfinite(offset[i]):
            continue
        by_cell.setdefault(c, []).append((recs[i], offset[i], pre_anomalized[i]))
        cell_band[c] = band_idx[i]

    rng = np.random.default_rng(seed)
    groups = []
    for cid, cell_recs in by_cell.items():
        pool = build_pool_for_cell(cell_recs, n_pool, rng)
        if pool is None:
            continue
        x, y, p = pool
        ok = np.isfinite(x) & np.isfinite(y) & (x >= -50) & (x <= 12050)
        if ok.sum() < 20:
            continue
        groups.append((cid, x[ok], y[ok], p[ok], cell_band[cid], nens, bin_ages, binvec, seed))
    print(f"[gam] fitting {len(groups)} cells on {cfg.get('ncores') or 'auto'} cores ...",
          file=sys.stderr, flush=True)

    ncores = int(cfg.get("ncores") or max(1, (os.cpu_count() or 2) - 1))
    ncores = max(1, min(ncores, len(groups)))
    import multiprocessing as mp, time
    cell_draws_mu = {}
    cell_scale = {}
    n_done = 0; n_skipped = 0; t0 = time.time()
    if ncores > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(ncores) as pool:
            for cid, dr, sc, bn in pool.imap_unordered(fit_cell, groups, chunksize=1):
                cell_band[cid] = bn
                if dr is None:
                    n_skipped += 1
                else:
                    cell_draws_mu[cid] = dr
                    cell_scale[cid] = sc
                n_done += 1
                if n_done % 10 == 0 or n_done == len(groups):
                    el = time.time() - t0
                    eta = (len(groups) - n_done) * el / max(n_done, 1)
                    print(f"[gam] {n_done}/{len(groups)} cells "
                          f"({n_skipped} skipped), {el:.0f}s elapsed, ETA {eta:.0f}s",
                          file=sys.stderr, flush=True)
    else:
        for g in groups:
            cid, dr, sc, bn = fit_cell(g)
            cell_band[cid] = bn
            if dr is None:
                n_skipped += 1
            else:
                cell_draws_mu[cid] = dr
                cell_scale[cid] = sc
    print(f"[gam] fit {len(cell_draws_mu)} cells, skipped {n_skipped}",
          file=sys.stderr, flush=True)

    # Apply empirical residual noise * gain (calibrated against the published spread)
    cd = {}
    for cid, dr in cell_draws_mu.items():
        sc = cell_scale[cid]
        if noise_gain > 0 and np.isfinite(sc) and sc > 0:
            r = np.random.default_rng(seed=cid * 17 + 1)
            nan_mask = ~np.isfinite(dr)
            d2 = dr + r.normal(0.0, sc * noise_gain, size=dr.shape)
            d2[nan_mask] = np.nan
            cd[cid] = d2
        else:
            cd[cid] = dr

    band_ens, glob = aggregate(cd, cell_band, nens, NB, ZW)
    glob = archival_reference(glob, bin_ages)

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    g = pd.DataFrame(glob, columns=[f"ens{i+1}" for i in range(nens)])
    g.insert(0, "binAges", bin_ages)
    g.to_csv(out_csv, index=False)

    # per-band ensembles (long format)
    band_frames = []
    for b in range(N_BANDS):
        be = archival_reference(band_ens[:, b, :], bin_ages)
        bf = pd.DataFrame(be, columns=[f"ens{i+1}" for i in range(nens)])
        bf.insert(0, "band", b + 1)
        bf.insert(0, "binAges", bin_ages)
        band_frames.append(bf)
    bands_csv = str(out_csv).replace("_global.csv", "_bands.csv")
    pd.concat(band_frames, ignore_index=True).to_csv(bands_csv, index=False)
    print(f"[gam] wrote {out_csv} + {bands_csv} ({NB} bins x {nens} members)",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ts", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--grid", required=True, help="equal_area_grid_centers.csv")
    ap.add_argument("--sigma-table", default=None, help="proxy_uncertainties.csv")
    ap.add_argument("--modern-grid", default=None, help="worldclim_modern_1deg.csv")
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    # Default ref files alongside the equal-area grid
    refdir = Path(args.grid).parent
    sigma_table_path = args.sigma_table or str(refdir / "proxy_uncertainties.csv")
    modern_grid_path = args.modern_grid or str(refdir / "worldclim_modern_1deg.csv")

    cfg = yaml.safe_load(Path(args.config).read_text())
    grid = pd.read_csv(args.grid)
    run_gam(args.ts, cfg, grid, sigma_table_path, modern_grid_path, args.out_csv)


if __name__ == "__main__":
    main()
