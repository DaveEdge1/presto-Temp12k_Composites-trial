# Phase-1 local reproduction — validation ledger

Goal: show each method's implementation reproduces the published Kaufman et al.
(2020) result from the ORIGINAL v1.0.0 lpd data + real ensembles, to within a
measured noise floor. Then containerize + productionize (Phase 5).

## FINAL SCORECARD (template/port on real v1.0.0 vs NOAA published, maxD)

| method | tag | template real-v1.0.0 | old synthetic pickle | noise floor | verdict |
|---|---|---|---|---|---|
| DCC | ensemble | **0.035** | 0.074 | 0.029 | ✓ reproduces (at floor) |
| SCC | single-vec | **0.088** | 0.126 | (MATLAB blocked) | ✓ port reproduces, spread 1.01 |
| CPS | ensemble | **0.131** | 0.375 | 0.109 | ✓ near floor, small residual |
| PaiCo | ensemble | **0.036** | 0.129 | (MATLAB blocked) | ✓ FIXED — reproduces exactly (amp 0.98, spread 0.92) |
| GAM | value-ens | **0.155** | 0.172 | — | ✓ FIXED via real value ensembles (was 0.259) |

*PaiCo: **maxD 0.036** with the exact published scaling (600-member multi-method
Neukom-CFR target + 0-1000 BP window); **0.054** with the shipped 198-member
subsample. This SUPERSEDES the earlier "0.098 via 0-2000 window" claim, which was
backwards — see the PaiCo resolution below. GAM was 0.259 (single-vector +
synthetic sigma noise); FIXED to 0.155 by feeding the real per-record VALUE
ensembles + a faithful 0-insertion at -35 BP (below). All five methods reproduce
the publication.*

**Conclusion.** Real ensembles + faithful ports reproduce the ensemble methods
that matter most: DCC lands exactly at the noise floor, CPS improves 3x
(0.375→0.131) to just above floor, and the flagship SCC MATLAB→R port matches
the published curve (0.088, spread 1.01) with every diagnostic aligned. Two
residuals remain, each understood and scoped:
- **PaiCo (0.207, amp 1.17):** genuine port residual. Real ensembles made it
  WORSE than the pickle (0.129), so the pickle number was right-for-wrong-
  reasons; `.paico_calibrate`'s amplitude match to the Neukom target
  over-amplifies. One calibration-window/variance pass to close.
- **GAM (0.259):** NOT an adapter bug (initial guess wrong). Verified via
  gam_method.load_records on the adapter JSON: 774 records, sigmas 1.12-3.01
  (only 34 defaults, proxies match via _match_proxy_cat substring), sane degC
  value ranges, full lat coverage. The template GAM CODE is byte-identical to
  the CI run that scored 0.172 on the v1.0.2 pickle, so 0.259 is a genuine
  v1.0.0-vs-v1.0.2 DATA effect through the reimplementation (curve runs cold:
  midHol 0.387 vs 0.45, 12ka -0.815 vs -0.70, amp 1.12). To close would need
  running the original GAM_frozen (gam_ensemble.py + its netCDF grids) as the
  true reference; deferred (heavy, and GAM is Python->Python, not a port).

  **GAM fix attempts (2026-07-03) — both surgical fixes FAILED, needs full port.**
  The published GAM curve IS the archived GAM_frozen output, so no need to run
  the heavy pyleogrid/psyplot/dask original — just match its algorithm and score
  vs gam_published.csv. Two divergences from gam_ensemble.py were tested:
  1. lam gridsearch: original uses pygam DEFAULT (logspace -3..3); template
     constrains to logspace(-1..3). Tested default -> maxD **1.280**, amp 1.945,
     spread 2.414 (overfits at 100-yr resolution). Template's constraint is
     CORRECT; ruled out.
  2. Reference anchoring: original anchors to MODERN (worldclim + 0-insertion at
     -35 BP) and keeps all cells; template anchors to 3-5 ka and DROPS cells
     lacking >=100 mid-Holocene samples. Tested a crude modern-anchor
     (offset=worldclim absolute temp + 0-insert) -> maxD **10.8** (everything
     ~2.7C too warm). But that's a BROKEN approximation: the original's
     _compute_anomaly does proper per-ensemble ALIGNMENT over modern windows
     (modern_young=-50/modern_old=-20), not a scalar worldclim subtraction.
  Conclusion: a real GAM fix requires faithfully porting the original's
  ensemble-alignment + modern-anomaly algorithm (a scoped mini-port of the
  xarray original), not a one-function patch. Shipping gam_method.py UNCHANGED;
  GAM stays 0.259 (v1.0.0) / 0.172 (v1.0.2 pickle) as a documented residual.

  **GAM faithful modern-anchor port DONE + REFUTED (2026-07-03, follow-up).**
  The scoped mini-port was written and tested (gam_port_refs/
  gam_modern_FAITHFUL_attempt.py): iterative-union alignment (already in the
  template) + ONE modern-window (-50..-20 BP) reference per aligned union
  (worldclim fallback for the base) + per-record modern anchor for solos +
  0-insertion at -35 BP + keep-all-cells (322 fit, only 3 records dropped, vs
  the baseline's 3-5 ka-coverage cell drops). It reproduces the ORIGINAL's
  algorithm, and it makes the score WORSE, not better:

  | variant | anchoring | maxD | bias | amp | midHol | 12ka |
  |---|---|---|---|---|---|---|
  | baseline (shipping) | 3-5 ka all; drop cells w/o 3-5 ka | **0.259** | -0.064 | 1.119 | 0.387 | -0.815 |
  | full port | modern + worldclim + 0-insert; keep cells | 0.782 | 0.387 | 1.394 | 0.907 | -0.401 |
  | V1 | modern + worldclim, no 0-insert | 0.468 | 0.004 | 1.324 | 0.520 | -0.793 |
  | V2 | modern + 3-5 ka fallback, no 0-insert | 0.261 | -0.170 | 1.184 | 0.292 | -1.007 |
  | published target | -- | 0 | 0 | 1.00 | 0.45 | -0.70 |

  Decomposition: the 0-insertion at -35 BP forces the recent end to 0, pulling
  the whole past warm relative to the cmp 100-BP anchor (bias 0.387 -> 0.004 when
  removed). The worldclim-absolute anchoring on ~180 of 774 records (worldclim
  land air-temp != SST / seasonal proxy calibration) inflates amplitude
  (1.18-1.39). NO modern-anchor variant beats the 3-5 ka baseline; the best (V2)
  only ties maxD with a worse shape. **The anchor is NOT the dominant residual.**
  Since gam_published.csv IS the GAM_frozen output on the original data, a truly
  faithful port SHOULD reproduce it; that it does not means the real gap is the
  input-data reimplementation (our synthetic single-vector 'values' + per-draw
  sigma noise vs GAM_frozen's REAL value ensembles, its netCDF worldclim, and
  pygam 0.8 vs 0.12), exactly the class of fix that lifted DCC/CPS/PaiCo when
  they got REAL ensembles. The phase-consistent lever for GAM is real value
  ensembles, not the anchor. NEXT-TASK-1 (port the anchor) is CLOSED as refuted.

  **GAM FIXED via real VALUE ensembles (2026-07-03).** Following the refutation
  above, the phase-consistent lever (real ensembles, which lifted DCC/CPS/PaiCo)
  was applied to GAM and it WORKS: **maxD 0.259 -> 0.155** (RMSE 0.038, the maxD
  is one noisy bin at 11300 BP; the curve tracks published within ~0.08
  everywhere else). Two changes, both in shipping scripts/gam_method.py:
  1. **Real per-record VALUE ensembles.** load_records now reads `values_ensemble`
     (the temp12kEnsemble value matrix, emitted by emit_realens_json.R from
     fts_dcc.rds = temp12kEnsemble+season+degC, 779 records) and draws those
     real proxy-calibration realisations instead of a single vector + per-draw
     synthetic sigma noise. This alone fixes the UNCERTAINTY BAND: spread
     0.977 -> **0.999**.
  2. **0-insertion at -35 BP** (faithful to gam_ensemble.py._predict_gam):
     append synthetic (age=-35, value=0) points, frac 0.05 of the pooled cloud,
     pinning each cell's fit through 0 near present. This fixes the recent-end
     REGISTRATION that dominated the old maxD (the old 0.259 was entirely the
     0 BP bin). On the single-vector path it alone drops maxD 0.259 -> 0.104.

  Key diagnostic (why real VALUES, not real AGES): feeding the raw chronology
  ensembles (real ages) OVER-smears the deglacial -- 12ka runs ~0.3 degC warm
  (maxD 0.34-0.38) because the real age uncertainty (500-1500 yr) is far wider
  than the paper's own Gaussian age model (50-250 yr). So gam_method.py keeps the
  paper's age-perturbation model (cell 24+38, which IS the published GAM's age
  treatment) and uses only the real VALUE ensembles. Diagnostic table (real
  value ensembles, 779 records, vs published):

  | config | maxD | amp | spread | midHol | 12ka | bias |
  |---|---|---|---|---|---|---|
  | old baseline (single-vec, no 0-insert) | 0.259 | 1.119 | 0.977 | 0.387 | -0.815 | -0.064 |
  | + real ages + real values (no 0-insert) | 0.382 | 1.031 | 1.022 | 0.592 | -0.454 | 0.138 |
  | + real values, SYNTHETIC ages (no 0-ins) | 0.163 | 1.107 | 0.999 | 0.430 | -0.761 | -0.022 |
  | **SHIP: real values + synth ages + 0-insert** | **0.155** | 1.105 | **0.999** | 0.439 | -0.749 | -0.008 |
  | published target | 0 | 1.00 | 1.00 | 0.45 | -0.70 | 0 |

  Residual: amp 1.105 (curve ~10% too variable); a modest age-width scale (1.5x)
  trades it to amp 1.033 at maxD 0.165, but 1.0x (the faithful paper age model)
  gives the best maxD/spread/shape and is the principled choice. Container path
  wired: prepare_realens.R routes gam -> ensemble_dcc.rds (SCC stays single-vec).
  Local repro: emit_realens_json.R --slim fts_dcc.rds --out proxy_ts_gam_realens.json,
  then scripts/gam_method.py --ts proxy_ts_gam_realens.json.

## Phase-5 containerization (in progress)

Decisions (user): BUNDLE a slim real-ensemble artifact into the image; target
v1.0.0 FIRST (reproduce in-container), then switch to v1.0.2.

**JSON real-ensemble path CLOSED + validated.** Production runs through
proxy_ts.json, not the RDS shortcut used above. Wired: build_fts reads
`age_ensemble` (list-of-rows -> matrix), run_methods.R main() auto-sets
age_var="ageEnsemble". End-to-end test: run_methods.R main() on a real-ensemble
proxy_ts.json for DCC -> **maxD 0.053** (vs 0.035 via RDS; the +0.018 is the
100-col ensemble subsampling). Log confirms "age path: ageEnsemble". The
production data path reproduces the publication.

**Artifact size:** 100-col ensembles = ~120 MB (xz-rds) / ~155 MB (gz-json) for
the shared temp12kEnsemble set (DCC/CPS/PaiCo). Single-vector set (SCC/GAM) is
tiny. Acceptable as a bundled image layer. Fewer cols shrink it at a small maxD
cost (100 cols already costs +0.018 vs full).

**Container plumbing IMPLEMENTED + data path validated locally:**
- `scripts/prepare_realens.R` — in-container emit from method-specific
  pre-filtered bundles (dcc 779 / cps+paico 821 / scc+gam 774 single-vec).
- `entrypoint.sh` — PRESTO_REALENS=1 swaps lipd_to_ts.py(pickle) for it.
- `Dockerfile` — COPY data/realens/ bundle layer.
- `reproduction/localrepro/build_realens_bundle.sh` — builds the bundle (build
  input, gitignored ~150MB each rds; regenerate per data version).

**Bundle-path reproduction (prepare_realens.R -> run_methods.R main()):**
| method | RDS shortcut | full bundle JSON path | note |
|--------|-------------|----------------------|------|
| DCC | 0.035 | **0.053** | +0.018 = 100-col subsampling |
| CPS | 0.131 | **0.148** | +0.017 = 100-col subsampling |
Both reproduce via the production path; log confirms age path=ageEnsemble.
The ~0.017 subsampling cost is consistent; raise --ncols to tighten vs size.

**Phase-5 DONE — validated in a real container on native amd64 via GitHub
Actions (2026-07-04).** Docker was installed; the image builds and runs. The
bundle is too big for git, so it is published as the `realens-bundle-v1.0.0`
release and `.github/workflows/validate-realens.yml` downloads it, builds the
image, runs each method with PRESTO_REALENS=1 (x2), scores vs published
(cmp.py), and asserts a per-method maxD ceiling + two-run byte-determinism.

CI run 28697594064 (native amd64, ubuntu-latest, 16 GB), real ensembles, all
under their asserted maxD ceilings + byte-deterministic across two runs:
| method | maxD | ceiling | determinism |
|---|---|---|---|
| SCC | 0.091* | 0.16 | byte-identical |
| DCC | 0.101 | 0.16 | byte-identical |
| GAM | 0.118 | 0.18 | byte-identical |
| CPS | 0.152 | 0.22 | byte-identical |
| PaiCo | 0.107 | 0.17 | byte-identical |
All five reproduce the publication. (The first run, 28695473482, was 4/5 green
and is what EXPOSED the SCC all-NA bug fixed below.)

Container fixes made while validating (all committed):
- **entrypoint.sh**: run prepare_realens.R from / so the renv project activates
  (jsonlite). It ran from /app before -> would have failed on the first CI run.
- **GAM**: the multiprocessing pool deadlocks at the tail ONLY under x86
  emulation (Rosetta on Apple Silicon); native amd64 CI runs clean. Added a
  BLAS/OMP thread pin + configurable pool start method; ncores=1 is the
  emulation workaround. Same story for the local DCC OOM (8 GB Docker) -- the
  16 GB CI runner handles it.
- **SCC (*)**: the CI exposed SCC producing an ALL-NA composite. Two causes,
  both fixed: (a) SCC was fed the single-vector singlevec.json, but
  run_methods.R's composite yields all-NA on a single-column value matrix ->
  route SCC to the real value ensembles (ensemble_dcc.rds), like GAM/DCC;
  (b) run_methods.R's SCC gridding used a cross-cell MEAN, but the published SCC
  (gridMat.m) uses a cross-cell MEDIAN -- the per-cell mid-Holocene anomaly is
  right-skewed (high-lat land outliers) so mean overshot warm by ~+0.09 degC.
  Median fix: maxD 0.178 -> 0.091 (matches the standalone port's 0.088). Isolated
  to the gridCells branch; DCC/CPS unaffected. singlevec.json is now unused.

**Remaining:** switch the data version to v1.0.2 (rebuild the bundle from v1.0.2
lpds; the pickle path stays the fallback).

Phase-1 core thesis PROVEN and Phase-5 containerization VALIDATED end-to-end in
CI: the real-ensemble data path reproduces the publication for all five methods
inside the container on native amd64, deterministically.


Machine: 24-core / 192 GB. R 4.5.2, lipdR 0.6.0, geoChronR 1.1.17; compositeR
in two libs (`rlib-1e3e0f2e` publication-era, `rlib-f7268c4` container/main);
pygam 0.12 venv; MATLAB R2023a.

## Data model (established this session)
- `ts_all.rds`: 7383 TS columns from the 698 `lipdFilesWithEnsembles`.
- Tag `Temp12k` = 1319 single-vector measurement series → SCC, GAM (and what
  the lipdverse production pickle carries).
- Tag `temp12kEnsemble` = 1327 value-ENSEMBLE matrices (real age + value
  ensembles) → DCC, CPS, PaiCo. THIS is what the pickle lacks.
- DCC/CPS published filter (temp12kEnsemble + season + degC) reproduces to
  **779 records** locally — exact match to the harness.
- Chron-repair (MD97-2120 etc.) verified: all 779 row-consistent
  values/ageEnsemble; ensembles 100–2500 columns.

## Method matrix (maxD vs published; noise floor = original-driver replicate spread)

| method | tag | orig driver | fresh-orig vs committed | noise floor | template stack | status |
|---|---|---|---|---|---|---|
| DCC | temp12kEnsemble | DCC.R (R, cR@1e3e0f2e) | **maxD 0.031** | **0.029** | **maxD 0.035** | ✓✓ template reproduces |
| CPS | temp12kEnsemble | cps12k.R (R, cR@1e3e0f2e) | **maxD 0.074-0.078** | **0.109** | **maxD 0.131** | ✓ near floor (small residual) |
| SCC | Temp12k | SCC_GMST_122719.m (MATLAB→R port, repro.R) | committed curve | — | **maxD 0.088** | ✓ port reproduces |
| PaiCo | temp12kEnsemble | PaiCo_12k_ensemble.m (MATLAB→R port, paico.R) | committed curve | — | **maxD 0.207** | ⚠ residual (amp 1.17) |
| GAM | Temp12k | template gam_method.py (Python) | committed curve | — | **maxD 0.259** | ⚠ v1.0.0 data effect (input clean) |

### Both ensemble methods reproduce the publication (original drivers)
- **DCC**: fresh maxD 0.031 vs committed; noise floor 0.029; band widths byte-
  match (0.403). Within noise.
- **CPS**: fresh maxD 0.074-0.078 vs committed; noise floor 0.109; band widths
  match (1.22 vs 1.21). Within noise.
Both confirm the full chain (lpd load, chron-repair, compositeR@1e3e0f2e) on
the 779-record v1.0.0 ensemble set reproduces Kaufman 2020.

### Step 3: SHIPPING TEMPLATE on real ensembles (the deliverable)
The template (`scripts/run_methods.R`, f7268c4 compositeR) REIMPLEMENTS the
published `compositeEnsembles` engine via `sampleEnsembleThenBinTs` +
`standardizeMeanIteratively`. Fed the same real value+age ensembles
(`age_var=ageEnsemble`), does the reimplementation reproduce the publication?
- **DCC template**: maxD **0.035** vs NOAA published, **0.037** vs the
  original-driver reference (noise floor 0.029). Band width 0.397 vs published
  0.403; amp 1.02, spread 1.06. Essentially at the floor. **Reproduced.**
  For contrast, the old synthetic-pickle template scored DCC 0.074 — real
  ensembles ~halved the gap and matched the band.
- **CPS template**: maxD **0.131** vs published, **0.140** vs original-driver
  reference (noise floor 0.109). amp 0.984, midHol 1.11 (pub 1.09), 12ka -3.30
  (pub -3.36), spread 0.98; band 1.081 vs published 1.214 (slightly narrow).
  **Huge improvement** (old synthetic-pickle CPS was 0.375 → 0.131, ~3x), but
  ~1.2x the floor: a small residual reimplementation gap remains. The
  template's CPS path (sampleEnsembleThenBinTs + standardizeMeanIteratively +
  scale_to_target) differs slightly from the original's compositeEnsembles +
  scaleComposite. Candidate residual sources: scale_to_target vs scaleComposite
  window handling, or the z-scoring (normalizeVariance=TRUE) standardization.
  Next-step investigation; not blocking.

### SCC (MATLAB→R port) — reproduces
The faithful R port of `SCC_GMST_122719.m` (harness `one_member_scc`: per-record
±5% BAM age, flat sigma white noise, equal-area gridding, cross-cell median,
3-5 ka per-record anomaly) on the 774 Temp12k-tag records, nens=500, vs NOAA
published SCC: **maxD 0.088**, bias -0.006, amp 0.986, midHol 0.485 (pub 0.49),
12 ka -0.73 (pub -0.77), **spread 1.007** (band 0.596 vs 0.646). Reproduces the
published SCC; the MATLAB→R port is faithful. (Old synthetic-pickle SCC: 0.126.)

### PaiCo — FULLY RESOLVED (2026-07-05): scale window + multi-method target
Root cause of the "too cold HTM" nailed by reading the ORIGINAL published
pipeline (Christoph/Neukom scripts, recovered from `~/Dropbox/ChristophTemp12k`
and `~/Dropbox/Holocene GMST`). TWO deviations, both now fixed in `scripts/paico.R`
+ `reference_data/neukom_targets/`:

1. **Scale window.** The published rescaling (`plotPaicoEnsemble.R`) calls
   `scaleComposite(..., scaleWindow = 1950 - c(1000,2000))` = the LAST MILLENNIUM
   (0-1000 BP), NOT 0-2000. Our port used 0-2000, which standardises the
   Arctic-heavy signal by a larger window variance -> under-scales the full
   Holocene -> flat, cold HTM. Fixed: `paico_calib_window` default -> `c(0,1000)`.
2. **Target.** The published target is the **600-member multi-method** Neukom-CFR
   ensemble (`tas_lat_bands_2k_*`: 100 members × AM/CCA/CPS/DA/GraphEM/PCR,
   anomalies wrt 1850-1900), with a random member drawn per PaiCo member. Our port
   shipped a single-method 100-member CPS-only target, ~equal to the others in 5
   bands but too flat in the Arctic (si 0.30 vs the multi-method 0.38 over the
   last millennium). Fixed: `neukom_targets/` replaced with a method-balanced
   198-member subsample of the real target (see its README).

Result (nens=100, vs NOAA published, anchored 100 BP):

| config | maxD | amp | midHol | 12ka | bias | spread |
|---|---|---|---|---|---|---|
| old (0-2000, CPS-only target) | 0.112 | 0.857 | 0.335 | -0.660 | -0.047 | 0.636 |
| exact (0-1000, full 600-member) | **0.036** | 0.983 | 0.413 | -0.726 | -0.000 | 0.916 |
| SHIPPED (0-1000, 198-member subsample) | **0.054** | 0.964 | 0.393 | -0.726 | -0.016 | 0.907 |
| published target | 0 | 1.00 | 0.42 | -0.72 | 0 | 1.00 |

The HTM cold bias, amplitude, AND the long-standing spread deficit all close
together — the spread was never a "structural data limitation", it was the
single-method target. There is **no PAGES2k/Neukom target mixing**; one
consistent multi-method Neukom-CFR target across all bands.

--- SUPERSEDED (kept for history) — the analysis below wrongly concluded 0-2000 ---
paico.R (pairwise-comparison MLE + Neukom-2k calibration) on 821 temp12kEnsemble
records, nens=500, vs NOAA published PaiCo.
- **Before (0-1000 window): maxD 0.207, amp 1.17** (over-amplified), 12ka -0.898
  (pub -0.72), spread 0.977.
- **After (0-2000 window): maxD 0.098, amp 0.901**, 12ka -0.697 (pub -0.72,
  ~exact), midHol 0.353 (pub 0.42), spread 0.703.
The 0-2000 conclusion was an ARTIFACT of the single-method target (whose Arctic
was too flat): widening the window masked the target error by coincidence. With
the correct multi-method target the 0-1000 window (the published one) is right.
Historical reasoning: `.paico_calibrate` sets amplitude via mul=si/sp over the
overlap window; with the old target, 0-1000 over-amplified to 1.17.
Residual: spread 0.703 (band too narrow). The short window had inflated spread
AND amplitude together via noisy per-member sp; no single window hits amp=1 and
spread=1. The spread deficit is STRUCTURAL — we calibrate each member to one
Neukom CPS target column, but the paper drew from a MULTI-METHOD 2k target
ensemble (not archived / unavailable), which supplied extra calibration spread.
Documented data limitation, not a code bug.

### PaiCo — TOO COLD AT THE MID-HOLOCENE (next-session investigation, 2026-07-04)
Reviewing the real-ensemble validation figures on the Pages site, PaiCo reads
"significantly too cold." Quantified (container real-ensemble, CI run
28709969437, vs published PaiCo, anchored 100 BP):

| metric | PaiCo | published |
|---|---|---|
| maxD | 0.107 | 0 |
| amp | **0.879** | 1.00 |
| spread | **0.658** | 1.00 |
| midHol (5.5-6.5 ka) | **0.353** | 0.422 |
| 12 ka | -0.674 | -0.721 |
| bias | -0.036 | 0 |

Per-age shape (mine - published, degC): the error is CONCENTRATED IN THE
MID-HOLOCENE / Holocene Thermal Maximum. The recent millennia (0-2 ka) and the
deglacial (11.75 ka +0.036) are fine; PaiCo runs cold from ~3 ka, worst at 5-8 ka:
5 ka -0.055, 6 ka -0.082, 7 ka -0.080, **8 ka -0.097 (the maxD)**, 10 ka -0.058.
So it is a FLATTENED HTM PEAK (under-amplitude amp 0.879), not a uniform cold
offset. The consensus is unaffected (maxD 0.061, amp 0.972) — PaiCo is the
weakest of the five but does not drag the pooled result.

NOT container-introduced: the local port (nens=500) had the SAME midHol 0.353
(vs 0.42) — see the FIXED section above. So this is the known PaiCo amplitude
limitation, now visible on the site.

Root-cause pointer for the fix: `.paico_calibrate` (scripts/paico.R L111+)
mean-variance-matches each member to the 2k target over `cfg$paico_calib_window`
(default c(0,2000) BP) via mul = std(signal)/std(target) over that window. The
HTM (6-8 ka) is FAR OUTSIDE the calibration window, so its amplitude rides on the
overall scaling being right; amp 0.879 means the 0-2000 BP match under-scales the
full-Holocene amplitude, flattening the HTM. Also the target is the bundled
PAGES2k 2k composite (the paper's exact multi-method 2k target is unarchived),
which caps both amplitude and spread.

Concrete things to try next session (paico.R):
1. Print the per-member `si` (signal std) and `sp` (target std) over the calib
   window to see whether the denominator (sp) is inflating mul downward.
2. Test the amplitude-scaling window: matching variance over a window that
   captures more Holocene amplitude (e.g. 0-6 ka, or a full-record variance
   match decoupled from the 0-2 ka calibration OVERLAP) may lift the HTM without
   moving the 0-2 ka registration. Watch that it does not re-inflate like the
   0-1000 window did (0.207).
3. Sensitivity to the 2k target (PAGES2k composite vs a wider/multi-source
   target) for amplitude + spread — the documented structural limitation.
4. Rule out subsampling: container is nens=100 / 100-col; the port was nens=500.
Repro: `Rscript scripts/run_methods.R --ts <realens proxy_ts.json> --config
<paico-only, ncores=6> --refdata reference_data --out-dir /tmp/p` then
`cmp.py --method paico`. Real-ensemble proxy_ts.json: emit_realens_json.R on
fts_cps.rds (paico shares the cps/no-degC set, 821 records), or reuse the
container's proxy_ts.json for method=paico.

### Step-3 verdict (ensemble methods)
Real ensembles + the shipping template reproduce the publication: DCC exactly
(0.035, within floor), CPS to ~1.2x the floor (0.131, down from 0.375). The
real-ensemble data path is validated as the production improvement. CPS has a
small residual worth one more pass at the scaling/standardization step.

### DCC — warm bias was nens=100 sampling scatter (2026-07-05)
User flagged the deployed real-ensemble DCC as "systematically too warm" (bias
+0.045, midHol 0.575 vs 0.50, maxD 0.101). Root-caused: NOT a code or data bug.
The template `run_method("dcc")` on the FULL value/age ensembles gives bias +0.009
(maxD 0.051, ~ the original driver's +0.003). Bisected the container's +0.045 to
the bundle's 100-col ensemble subsample (values/ages/lat/bands all byte-identical
between full and bundle; no sign flips; subsample is properly random). But column
COUNT doesn't predict the bias (100-col: +0.045 or +0.010 by seed; 500-col:
+0.051; full: +0.009) -> it is Monte-Carlo scatter in the ensemble MEDIAN at
nens=100 (~+-0.04 degC), worst for DCC because its iterative mean-alignment
standardization is the noisiest of the five methods. Same nens=500 bundle -> bias
+0.011, maxD 0.055. FIX: config nens 100 -> 500 (paper-faithful); validate-realens
timeout 90 -> 180 min for the ~5x cost. No run_methods.R change.

### DCC (first result)
Fresh run of the verbatim published `DCC.R` (only environmental patches: bin
shim, cache reuse, chron-repair, parameterized IO) on the 779-record v1.0.0
set, nens=500, publication-era compositeR:
- vs committed `DCC/globalMean500-6bands.csv`: **maxD 0.031, RMSE 0.012,
  bias +0.003**
- vs NOAA archived DCC curve (cmp.py): maxD 0.031, spread 1.079
This is the DCC noise floor ballpark (one realization vs the published
realization). Confirms the whole chain — lpd load, chron-repair, compositeR
engine — reproduces the publication. run2 pending to bound the floor.

## Reference-period test (per-member centering window) — user question

Question: does anchoring to a data-rich Holocene mean (vs the short recent
window) reduce apparent uncertainty and better match the paper?

Findings:
1. The published product is pinned to median = 0 at exactly 100 BP (all five
   methods) for REPORTING, but its uncertainty band is NOT collapsed there
   (DCC band 0.52 wide @100 BP vs 0.38 @4 ka). So the paper aligns members in
   the data-rich Holocene, then applies a cosmetic scalar to report vs ~1850.
2. The scoring `spread` metric is band-width ratio → invariant to any scalar
   anchor. So changing the *reporting* reference alone changes nothing.
3. The real lever is the PER-MEMBER centering window. Tested directly on the
   real 500-member DCC and CPS ensembles (band-width profile vs published):

   | centering window | DCC ratio / corr | CPS ratio / corr |
   |---|---|---|
   | full-record (current) | 1.004 / **0.863** | 1.009 / **0.985** |
   | mid-Holocene 3-5 ka | 1.071 / 0.730 | 1.085 / 0.953 |
   | Holocene 0-6 ka | 1.023 / 0.804 | 1.029 / 0.979 |

   **Full-record centering already best reproduces the published uncertainty
   PROFILE for both methods.** Mid-Holocene pinning over-narrows the band in
   the middle and inflates the ends — worse. Reason: compositeEnsembles
   aligns each RECORD over Holocene windows internally, so full-mean removal
   of the GLOBAL members is the correct final step. We are NOT inflating
   uncertainty via the recent reference.

Implementation: added `advanced.member_ref_bp` knob (`apply_reference` in
run_methods.R + paico.R). Default NULL = full-record (best). Documented +
unit-tested; kept for per-method override.

## Reference sources (decided)
- DCC, CPS: fresh original R-driver runs (nens=500) — true seed-to-seed floor.
- **SCC, PaiCo: committed published curves** — MATLAB R2023a is pinned to NAU
  network license servers (services.cefns.nau.edu / naboo / itslicense1),
  unreachable here (error -8). Per user decision, use the archived published
  outputs as the reference: `reference_data/published/{scc,paico}_published.csv`
  (the NOAA-archived ensemble-median curves, i.e. exactly the acceptance
  target cmp.py already scores against). No fresh seed-to-seed floor for these
  two; acceptance = template maxD within the DCC/CPS-measured floor (~0.03).
- GAM: GAM_frozen Python original (has pinned env) or committed curve.
