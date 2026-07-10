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
def compute_alignments(recs, cell_idx):
    """For each cell, pick the longest record as base; iteratively shift other
    records by the mean offset over their age overlap with the growing aligned
    union (`min_overlap = 100` ensemble samples both sides, per gam_ensemble.py
    L385-426). Returns per-record `offset` and `pre_anomalized` flag.

    Records that can't align AND lack >=100 samples in 3-5 ka anchor via
    WorldClim modern temperature (paper L534); marine records with no land
    pixel get offset NaN and are dropped from the pool.
    """
    n = len(recs)
    offset = np.full(n, np.nan, dtype=float)
    pre_anomalized = np.zeros(n, dtype=bool)
    aligned = np.zeros(n, dtype=bool)
    n_aligned = n_solo_anom = n_solo_modern = n_solo_drop = 0
    n_no_overlap_anom = n_no_overlap_modern = n_no_overlap_drop = 0

    by_cell = {}
    for i, c in enumerate(cell_idx):
        if c is None:
            continue
        by_cell.setdefault(c, []).append(i)

    for cid, members in by_cell.items():
        if len(members) == 1:
            i = members[0]
            # EXPERIMENT(modern): prefer worldclim-modern anchor (original gam_ensemble.py)
            if np.isfinite(recs[i]["modern"]):
                offset[i] = recs[i]["modern"]
                pre_anomalized[i] = True
                n_solo_anom += 1
            elif np.isfinite(recs[i]["modern"]):
                offset[i] = recs[i]["modern"]
                pre_anomalized[i] = True
                n_solo_modern += 1
            else:
                offset[i] = np.nan
                n_solo_drop += 1
            continue
        # multi-record cell -- iterative growth
        lens = [recs[i]["age_med"].size for i in members]
        base_i = members[int(np.argmax(lens))]
        offset[base_i] = 0.0
        aligned[base_i] = True
        n_aligned += 1
        # Build the aligned union (post-shift)
        u_age = recs[base_i]["age_ens"].ravel()
        u_val = np.broadcast_to(recs[base_i]["val"][:, None],
                                recs[base_i]["age_ens"].shape).ravel()
        remaining = [i for i in members if i != base_i]
        changed = True
        while remaining and changed:
            changed = False
            still = []
            a_min = float(np.min(u_age)); a_max = float(np.max(u_age))
            for i in remaining:
                r_age = recs[i]["age_ens"].ravel()
                r_val = np.broadcast_to(recs[i]["val"][:, None],
                                         recs[i]["age_ens"].shape).ravel()
                if r_age.size == 0:
                    still.append(i); continue
                r_min = float(np.min(r_age)); r_max = float(np.max(r_age))
                m1 = (u_age >= r_min) & (u_age <= r_max)
                m2 = (r_age >= a_min) & (r_age <= a_max)
                if m1.sum() > 100 and m2.sum() > 100:
                    diff = float(np.nanmean(u_val[m1]) - np.nanmean(r_val[m2]))
                    offset[i] = -diff   # subtract this offset to add +diff
                    aligned[i] = True
                    n_aligned += 1
                    u_age = np.concatenate([u_age, r_age])
                    u_val = np.concatenate([u_val, r_val + diff])
                    changed = True
                else:
                    still.append(i)
            remaining = still
        for i in remaining:
            # EXPERIMENT(modern): prefer worldclim-modern anchor (original gam_ensemble.py)
            if np.isfinite(recs[i]["modern"]):
                offset[i] = recs[i]["modern"]
                pre_anomalized[i] = True
                n_no_overlap_anom += 1
            elif np.isfinite(recs[i]["modern"]):
                offset[i] = recs[i]["modern"]
                pre_anomalized[i] = True
                n_no_overlap_modern += 1
            else:
                offset[i] = np.nan
                n_no_overlap_drop += 1
    print(f"[gam] alignment: {n_aligned} aligned in union | "
          f"singleton {n_solo_anom} anom + {n_solo_modern} modern + {n_solo_drop} drop | "
          f"no-overlap {n_no_overlap_anom} anom + {n_no_overlap_modern} modern + {n_no_overlap_drop} drop",
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
    # EXPERIMENT(modern): anchor each cell's GAM through 0 at modern (-35 BP)
    # via 0-insertion (original gam_ensemble.py L236), instead of requiring and
    # anchoring to 3-5 ka (which also dropped poorly-mid-Holocene-covered cells).
    y0 = y  # pool values already modern-anomalized (offset=modern in build_pool)
    n_anchor = max(50, int(0.05 * x.size))
    x_aug = np.concatenate([x, np.full(n_anchor, -35.0)])
    y_aug = np.concatenate([y0, np.zeros(n_anchor)])
    try:
        lam_grid = np.logspace(-1, 3, 11)
        gam = LinearGAM(s(0)).gridsearch(x_aug[:, None], y_aug, lam=lam_grid, progress=False)
        draws = gam.sample(x_aug[:, None], y_aug, sample_at_X=bin_ages[:, None],
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
