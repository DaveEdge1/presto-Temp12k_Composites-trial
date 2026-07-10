# Local Phase-1 reproduction — session handoff (START HERE)

Authoritative state for the LOCAL real-ensemble reproduction work. Supersedes
the CI-era `reproduction/ci/HANDOFF.md` (that phase's determinism/CPS-lever work
is DONE and merged to main; see it only for history). Full ledger:
`reproduction/localrepro/VALIDATION.md`. Memory: `project_phase1-localrepro.md`.

Branch: **`local/phase1-real-ensembles`** (all work pushed; not merged to main).

## What this phase proved
Feed the SHIPPING template real per-record age+value ensembles (from the 698
v1.0.0 lpd files) instead of the lipdverse pickle's synthetic single vectors,
and it reproduces the published Kaufman 2020 curves. Scorecard (template/port on
real v1.0.0 vs NOAA published, maxD; noise floors DCC 0.029 / CPS 0.109):

| method | maxD | status |
|--------|------|--------|
| DCC   | 0.035 | ✓ at floor |
| SCC   | 0.088 | ✓ port 0.088 / container run_methods.R 0.091 (cross-cell median) |
| PaiCo | 0.098 | ✓ MATLAB→R port; FIXED (was 0.207, calib-window bug) |
| CPS   | 0.131 | ✓ near floor |
| GAM   | 0.155 | ✓ FIXED via real value ensembles + 0-insert (was 0.259) |

Old synthetic-pickle scores for contrast: DCC 0.074, SCC 0.126, CPS 0.375,
PaiCo 0.129, GAM 0.172.

## Environment (this machine — 24 core)
- Temp12k clone: `~/GitHub/Temperature12k` (698 lpds in
  `ScientificDataAnalysis/lipdFilesWithEnsembles`; original drivers DCC.R,
  cps12k.R, SCC_GMST_122719.m, PaiCo/*, GAM_frozen/scripts/gam_ensemble.py).
- compositeR two libs: `reproduction/localrepro/rlib-1e3e0f2e` (publication),
  `rlib-f7268c4` (container/main). Select via `R_LIBS_USER=...`.
- pygam venv: `reproduction/localrepro/venv`.
- Caches (gitignored, rebuild if missing): `reproduction/localrepro/cache/`
  — `ts_all.rds` (7383 TS cols), slim caches `fts_{dcc,cps,scc}.rds`,
  `proxy_ts_gam.json`. Rebuild: `bash reproduction/localrepro/run_repro.sh cache`
  then `Rscript reproduction/localrepro/build_slim.R --method cps ...`.
- **MATLAB R2023a is license-BLOCKED** (NAU network servers, error -8) → SCC/
  PaiCo use the committed NOAA published curves as reference, not fresh reruns.

## Key scripts (reproduction/localrepro/)
- `run_repro.sh cache|harness <m>|score` — driver.
- `originals/orig_{dcc,cps}.R` — verbatim published R drivers (patched only for
  env: bin shim, cache reuse, chron-repair, IO). Reproduce refs within noise.
- `run_template_realens.R` — SHIPPING run_method on the slim cache (RDS path).
- `emit_realens_json.R` / `prepare_realens.R` — emit real-ensemble proxy_ts.json.
- `build_realens_bundle.sh` — builds the container bundle (data/realens/).
- `export_ts_mat.R` + `ndjson_to_tsmat.py` + `run_scc.m` — MATLAB SCC path (blocked).

## Shipping-code changes made this phase (on the branch)
- `scripts/run_methods.R`: `build_fts` reads `age_ensemble`; `run_method` takes
  `cfg$age_var` and main() auto-sets it to "ageEnsemble" when real ensembles are
  present. `apply_reference` gains `member_ref_bp` (default full-record = best).
  **SCC gridding uses cross-cell MEDIAN** (was mean; gridCells branch only, so
  DCC/CPS unaffected) — the SCC fix, maxD 0.178→0.091.
- `scripts/paico.R`: `cfg$paico_calib_window` default **c(0,2000)** — the PaiCo
  fix (0.207→0.098).
- `scripts/gam_method.py`: **GAM FIXED.** load_records reads `values_ensemble`
  (real value ens; sigma=0 when present); fit_cell adds 0-insertion at -35 BP
  (advanced.gam_zinsert_frac, default 0.05); BLAS/OMP thread-pin + configurable
  pool start method (advanced.gam_mp_context / GAM_MP_CONTEXT). Backward-compat:
  single-vector JSON still works (synthetic sigma).
- `scripts/prepare_realens.R`: **all 5 methods use real ensembles** — scc/dcc/gam
  → ensemble_dcc.rds, cps/paico → ensemble_cpspaico.rds. singlevec.json now unused.
- `entrypoint.sh`: `PRESTO_REALENS=1` mode; runs prepare_realens.R **from /** so
  renv/jsonlite activate (was broken: ran from /app).
- `Dockerfile`: COPY `data/realens/` bundle layer.
- `.github/workflows/validate-realens.yml`: per-method matrix, PRESTO_REALENS=1
  x2, score vs published + maxD-ceiling assert + byte-determinism.

## STATUS: Phase-1 (fidelity) + Phase-5 (containerization) BOTH DONE.
All 5 methods reproduce the published Kaufman 2020 curves from real v1.0.0
ensembles, validated in a real container on native amd64 CI (byte-deterministic).
GAM fix + Phase-5 details are in VALIDATION.md; the two headline results:
- **GAM**: maxD 0.259→0.155 (real VALUE ensembles + 0-insertion; the
  modern-anchor port was tried and REFUTED first). Test:
  `emit_realens_json.R --slim cache/fts_dcc.rds --out /tmp/g.json --ncols 100`
  then `gam_method.py --ts /tmp/g.json ...` → cmp.py maxD 0.155.
- **Phase-5 CI**: run 28697594064 green — SCC 0.091, DCC 0.101, GAM 0.118,
  CPS 0.152, PaiCo 0.107, all byte-identical. Bundle published as the
  `realens-bundle-v1.0.0` release (gh needs `--repo nickmckay/...`; pushing
  workflow files needs the token `workflow` scope).

## REMAINING WORK (for the next session)
1. **PaiCo too cold at the mid-Holocene (user-flagged 2026-07-04).** On the
   Pages validation figures PaiCo reads too cold: it FLATTENS the Holocene Thermal
   Maximum (5-8 ka ~0.08-0.10 degC too cold; maxD at 8 ka), amp 0.879, spread
   0.658. Recent millennia + deglacial are fine, so it's under-amplitude, not a
   cold offset (consensus is fine: maxD 0.061). Same as the local port residual
   (midHol 0.353), so it's the known PaiCo amplitude limit. Root cause + concrete
   experiments in VALIDATION.md "PaiCo — TOO COLD AT THE MID-HOLOCENE": start at
   `.paico_calibrate` (scripts/paico.R L111+), the mul=si/sp variance match over
   `paico_calib_window` c(0,2000) — the HTM is outside that window. **This is the
   next-session task the user wants.**
2. **v1.0.2 production data version.** Rebuild the bundle from v1.0.2 lpds
   (`build_realens_bundle.sh` after rebuilding the fts caches from v1.0.2), re-run
   the CI validation, confirm scores hold. Pickle path stays the fallback.
3. **Merge decision.** Branch `local/phase1-real-ensembles` is ready for review:
   guarded/backward-compatible shipping-script changes + the container
   real-ensemble path + CI workflows + all localrepro tooling. Nothing merged to
   main yet.
4. **Other residuals (low priority):** CPS 0.131 (small reimpl gap), GAM amp 1.105.

## GitHub Pages diagnostics site (added 2026-07-04)
Live at **http://nickmckay.org/presto-Temp12k_Composites-trial/** — a landing
page with a tile per diagnostic: pickle-path validation, **real-ensemble
validation (v1.0.0)**, and the spatial visualizer. Each validation report overlays
the reconstructed GMST (median + 5-95% band) vs published Kaufman 2020, per method
+ consensus. Wiring: `validate-realens.yml` has a `combine` job that runs the
container combine stage on the 5 real-ensemble methods and uploads
`realens-diagnostics-<run_id>`; `visualize.yml` (the single Pages deployer) fetches
the latest such artifact into `/realens/` and rebuilds a unified landing page.
To REFRESH the real-ensemble tile after a new validate-realens run: dispatch
`visualize.yml` (`gh workflow run visualize.yml --ref <branch>`); it grabs the
latest artifact (not auto-triggered by validate-realens — decoupled on purpose).
NOTE: Pages deploys are branch-gated to the github-pages environment policy; the
feature branch was added to the allowed list to deploy from the branch (remove it
or rely on merge-to-main later). Pages was enabled with source=GitHub Actions.
