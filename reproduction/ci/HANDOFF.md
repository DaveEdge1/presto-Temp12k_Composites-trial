# Fidelity work — session handoff (for continuing on another machine)

> **SUPERSEDED (2026-07-03).** This CI-era phase (determinism fixes, CPS lever
> sweep) is DONE and MERGED to `main`. Active work moved to LOCAL real-ensemble
> reproduction — **start at `reproduction/localrepro/HANDOFF.md`** (branch
> `local/phase1-real-ensembles`). Keep this file for CI-loop mechanics/history.

Goal: improve the presto-Temp12k container's fidelity to Kaufman et al. (2020),
verified via **GitHub Actions CI** (the container + lipdverse pickle + methods run in
the cloud). No local Docker/R/pandas required to drive it.

## New-machine requirements (lightweight driver only)
- `git`; `gh` CLI **authenticated** (`gh auth login`; needs `repo`+`workflow` scopes on
  the fork `nickmckay/presto-Temp12k_Composites-trial`).
- `python3` + `numpy` (scoring only — NO Docker, R, or pandas needed here).
- Claude Code. Start a fresh session in the cloned repo; have it read this file + the
  project memory (see below).
- **One-time:** Actions must be enabled on the fork (Actions tab → "enable workflows").
  Already done on the current fork.

## How the loop works
1. Make a single-variable change on a branch off `main`; push to `origin`.
2. `bash reproduction/ci/dispatch_all.sh <branch>` — dispatches `reconstruct.yml`
   (archived payload Temp12k v1_0_2). CI runs 5 methods + `combine`, which **commits
   `results/methods/*.csv` back to the branch**.
3. `bash reproduction/ci/score_branch.sh <branch>` — fetches + scores each method's maxD/
   bias/spread. `cmp.py` reproduces the authoritative `comparison.json` EXACTLY.
   Each CI run ≈ 25-30 min; runs on different branches go in parallel.

## Current results (maxD vs published; lower better; nens=100 noise floor ≈ ±0.05)
| method | committed baseline (STALE) | fresh `main` | `exp/scc-median` | `exp/cr-pin-1e3e0f2e` |
|---|---|---|---|---|
| SCC | 0.225 | 0.144 | 0.184 | 0.127 |
| DCC | 0.200 | 0.163 | 0.139 | 0.085 |
| GAM | 0.167 | 0.165 | 0.168 | 0.169 |
| CPS | 0.656 | 0.422 | 0.370 | 0.276 |
| PaiCo | 0.115 | 0.140 | 0.096 | 0.167 |

### Conclusions
1. **Committed `comparison.json` is STALE** — fresh `main` is CPS 0.422 (not 0.656),
   better across the board (base image `lipdbase2` is a mutable tag / pickle drifted
   favorably). Regenerate the repo's headline baseline from a fresh `main` run.
2. **compositeR publication pin (1e3e0f2e) is a real win** vs fresh main: CPS −0.15,
   DCC −0.08; GAM flat (Python control) and PaiCo within noise (its core is the recent
   R port, uses compositeR only for binning). Confirm at nens=500, watch PaiCo.
3. **`exp/scc-median` is NOT a win** (SCC 0.144→0.184, within noise; gridMat.m absent so
   unverifiable). Drop it.
4. **CPS record-subset was a dead end** — template keeps 809 records vs harness 779
   (~30 apart), and the pickle has no `temp12kEnsemble` tag (uses `inCompilation="Temp12k"`).
   The CPS gap is the **age-ensemble + value-ensemble** axes, not subset. ~~The template
   propagates NO age uncertainty (`ageVar="age"` single vector; README's "BAM ±5%" claim
   is false)~~ **CORRECTED 2026-07-02: WRONG — compositeR's sampleEnsembleThenBinTs
   auto-runs BAM (±5%) on single age vectors; age uncertainty IS propagated and
   README is correct. See "Next steps" item 3 below.**
5. **Noise floor ≈ ±0.05 at nens=100.** Small effects need nens=500 or replicates.

## DETERMINISM ACHIEVED (2026-07-02, branch `exp/seed-rng`) — merge candidate

Back-to-back CI runs 28619537951 / 28621091018 on `exp/seed-rng` produced
byte-identical results/ (all five method CSVs, reconstruction, comparison.json;
run 6's auto-commit had nothing to commit). Three fixes, all on that branch:
1. `scripts/run_methods.R` + `scripts/paico.R`: per-member, per-method
   `set.seed` from `advanced.seed` (default 42) — R methods had NO seeding.
2. `scripts/lipd_to_ts.py`: value-ensemble RNG was seeded with builtin
   `hash(rid)`, which is SALTED PER PROCESS (PYTHONHASHSEED) — every previous
   run fed different value ensembles to ALL methods. Now `zlib.crc32(rid)`.
   This was the dominant noise source.
3. `scripts/gam_method.py`: pygam's `gam.sample()` uses numpy's legacy GLOBAL
   RNG in forked workers, unseeded — now seeded per cell.

**Canonical seeded scores (seed=42, nens=100, main + fixes), maxD vs published:**
| SCC | DCC | GAM | CPS | PaiCo |
|---|---|---|---|---|
| 0.126 | 0.074 | 0.172 | 0.375 | 0.129 |

All future A/B comparisons against these are exact paired deltas (any nonzero
delta is caused by the change; judge whether it generalizes by re-running the
pair at a second seed if the delta is small).

## Next steps (priority) — revised 2026-07-02 evening
1. ~~Merge `exp/seed-rng` to main~~ **DONE 2026-07-02 ~21:40 UTC** (fast-forward
   to 4f96457, user-approved). Main now carries the seeded pipeline AND the
   verified seeded results; `results/validation/comparison.json` on main matches
   the canonical seeded baseline (stale CPS-0.656 headline replaced).
2. compositeR pin paired re-test: **DONE — PIN REJECTED (2026-07-02 ~22:35 UTC,
   run 28624661235 on `exp/cr-pin-seeded`).** Paired deltas vs canonical seeded
   baseline: SCC +0.024, DCC exactly 0, GAM exactly 0, CPS +0.075, PaiCo +0.044
   (positive = worse). The publication-era compositeR (1e3e0f2e) REGRESSES
   CPS/PaiCo/SCC vs the container's f7268c4; DCC's code path is identical
   between versions and GAM is the no-compositeR control (both exactly 0,
   which also confirms determinism holds across a different Docker build).
   The pin line is CLOSED. (Single-seed caveat: deltas are exact at seed 42.)
3. ~~CPS/SCC age-uncertainty propagation~~ **MOOT — AUDIT FINDING CORRECTED
   2026-07-02.** The 2026-07-01 audit claimed the template propagates no age
   uncertainty and the README "BAM ±5%" claim is false. WRONG: compositeR's
   `sampleEnsembleThenBinTs` (both container f7268c4 AND publication 1e3e0f2e)
   automatically runs `geoChronR::simulateBam` (bernoulli, param=0.05 = ±5%)
   per member whenever `ts[[ageVar]]` is a single vector — which is what
   build_fts passes. There is no inner tryCatch: had simulateBam errored, every
   record would bin to NA and results would be empty; they aren't. So BAM age
   uncertainty IS propagated per member per record, and README_NOTES.md is
   CORRECT as written. The CPS gap (seeded maxD 0.375) must come from another
   axis: value-ensemble regeneration quality (10-col AR1 from single vector vs
   paper's real ensembles), CPS scaling target, or record subset.

## CPS lever sweep (2026-07-02 ~23:40 UTC): ALL THREE REJECTED

Three single-variable paired experiments off seeded main (canonical baseline
SCC 0.126, DCC 0.074, GAM 0.172, CPS 0.375, PaiCo 0.129). In every run the
untouched methods scored EXACTLY baseline (perfect controls):

| experiment | branch / run | CPS | other deltas | verdict |
|---|---|---|---|---|
| A: VALUE_ENSEMBLE_SIZE 10→100 | `exp/vens100` / 28627452408 | 0.547 (+0.172) | SCC +0.066, DCC +0.015, PaiCo +0.001, GAM 0 | REJECTED |
| B: CPS scale to Neukom targets | `exp/cps-neukom-target` / 28627453339 | 3.168 (catastrophic; spread 0.17) | all others exactly 0 | REJECTED (target units/variance mismatch) |
| C: CPS degC-only records (faithful to cps12k.R) | `exp/cps-degc` / 28627454180 | 0.624 (+0.249) | all others exactly 0 | REJECTED |

**Readings:**
- Main's current config is locally optimal on all three axes; the CPS 0.375
  baseline stands.
- Counterintuitive pattern: every change MORE faithful to the published
  pipeline (publication compositeR pin, degC-only CPS subset) scores WORSE
  against the published curve. The template's divergences apparently
  compensate for each other (or for pickle-data drift). Single-axis
  faithfulness restoration is a dead approach; remaining gap likely needs
  either multi-axis simultaneous changes or is data-path (pickle vs
  publication input) and unfixable from here.
- A's dose-response (10 cols good, 100 cols bad) suggests trying
  VALUE_ENSEMBLE_SIZE=1 (or dropping regeneration so compositeR simulates
  per-member AR1 noise from paleoData_uncertainty1sd) as a cheap follow-up.

## Follow-up: exp/no-vens (run 28629630608, 2026-07-03 ~00:37 UTC) — REJECTED

Dropped regenerated ensembles entirely (compositeR simulates fresh per-member
AR noise from uncertainty1sd, restoring per-method ar: SCC 0, DCC/CPS sqrt(.5)).
Paired deltas: SCC +0.045, DCC +0.014, GAM exactly 0 (control), CPS +0.089,
PaiCo +0.061. WORSE across the board.

**Value-noise dose-response is now bracketed and non-monotonic:**
fresh-per-member (max diversity) CPS +0.089 | 10 pre-baked cols = 0.375 BEST |
100 pre-baked cols CPS +0.172. Reusing 10 noise realizations across 100
members effectively averages down the injected variance; both more noise
diversity AND the published per-method noise structure score worse. The
untested limit is ZERO added value noise (emit 2 identical base columns so
compositeR's NCOL>1 path always draws the clean vector) — would tell whether
value noise helps at all, though spread (0.86 SCC) would drop further.
VALUE_ENSEMBLE_SIZE=10 stays. The CPS 0.375 baseline stands.

## Data-path drift quantified (2026-07-03, local analysis) — DRIFT IS STRUCTURAL, NOT INVENTORY

Method: downloaded the CI proxy artifact (lipd_legacy.pkl, run 28619537951) and
diffed its record set against the publication's actual ensemble input list
(nickmckay/Temperature12k → ScientificDataAnalysis/lipdFilesWithEnsembles,
698 .lpd files = the datasets the published DCC/CPS/PaiCo consumed).

**Inventory drift: negligible (~1.3%).**
- Pickle Temp12k tag: 1319 records / 696 datasets. Case-insensitive overlap
  with the publication's 698: 696. (An earlier exact-case diff suggesting 47
  missing was filename-case noise.)
- Truly missing: 2 datasets — `Duranunlak.EPD` (absent), and
  `Gunnarsfjorden.Allen.2007` (PRESENT in the pickle but tagged lowercase
  `temp12k`, so the exact-match tag filter drops it — fixable hygiene bug:
  match inCompilation case-insensitively).
- 7 more datasets dropped because v1_0_2 reclassified their seasonalityGeneral
  to `summer+`/`winter+` (excluded by both published and template filters).
- Pickle-only datasets in the ensemble sense: 0.
- Selection funnel reproduced locally: 1319 → 807 season-ok (log: 809) → 766
  degC; 721/807 carry a stated temperature12kUncertainty.

**Structural drift: the real gap.** The publication's 698 lpd files exist to
carry REAL per-record age + value ensembles (Bacon-style chronologies); the
lipdverse pickle collapsed everything to single vectors. The template
approximates with BAM ±5% ages + 10-col AR1 values — and today's sweep showed
every perturbation of that approximation scores worse. Conclusion: the
residual CPS gap is dominated by real-vs-synthetic ensemble structure, which
CANNOT be recovered from the pickle. The one heavy-but-concrete path: have CI
prepare-data fetch the 698 publication lpd files and attach their real
ensembles (ageEnsemble matrices; compositeR's NCOL>1 path consumes them
natively). Large change; park unless per-method fidelity matters enough.

## GAM/CPS diagnostic round 1 (runs 2863155xxxx, 2026-07-03 ~01:30 UTC)

Paired deltas vs canonical baseline; untouched methods exactly baseline in
every run (controls clean):
- `exp/gam-lam0` (lam floor 0.1→1): GAM maxD 0.172→**0.165** (−0.007), amp
  unchanged (1.136), spread 0.978. **First WIN of the campaign** (small,
  exact). amp overshoot is lam-insensitive in this range. Dose step 2
  (floor→10, `exp/gam-lam10`) dispatched.
- `exp/gam-ng05` (noise_gain 1.0→0.5): GAM maxD 0.153 (−0.019) BUT spread
  collapses 0.977→0.489. noise_gain is the spread knob, already calibrated;
  REJECTED as a maxD lever (trade-off, not a win). Keep 1.0.
- `exp/cps-norescale` (rescale=FALSE): CPS **bit-identical to baseline** —
  a perfect no-op. Mechanism understood: apply_reference subtracts member
  means and re-anchors at 100 BP, so any mean-shift from rescale cancels.
  IMPORTANT: CPS's +0.2 warm offset is therefore a SHAPE error (late-Holocene
  decline into the anchor too shallow), NOT a level error. Mean-matching
  hypotheses are untestable/irrelevant under referencing.
- `exp/cps-scalewin1k` (scale window 0-1000 BP): CPS 1.076 — REJECTED
  decisively; full-2k window stays (knob re-deadened; consider removing it
  from config docs instead).

## Round 2 (runs 28636484065 / 28636484989, 2026-07-03 ~03:55 UTC)

- `exp/gam-lam10` (lam floor→10): GAM 0.189 (+0.017) — over-smoothed.
  **lam dose-response COMPLETE: floor 0.1 → 0.172 | floor 1 → 0.165 BEST |
  floor 10 → 0.189. PROMOTE `exp/gam-lam0` (floor=1) to main.**
- `exp/scc-flatsigma`: bit-identical no-op. Mechanism: compositeR's
  uncVar/defaultUnc noise path only fires when paleoData_values is
  single-column; the template's 10-col value ensembles bypass it entirely, so
  SCC's uncertainty model knobs are DEAD CODE on this pipeline. SCC spread
  (0.864) not tunable from here; would need per-method value handling. PARKED.

## Branch map (all pushed to origin)
- `overnight/fidelity-plan` — housekeeping commits (sin-lat weights, dead-knob docs),
  audit refresh, `reproduction/audit/overnight_2026-07-01.md` write-up, these ci helpers.
- `exp/cps-t12kensemble` — FAILED (temp12kEnsemble tag absent); dead end.
- `exp/scc-median` — SCC rowMeans→median; not a win.
- `exp/cr-pin-1e3e0f2e` — compositeR publication pin; PARKED (noise).
- `baseline/fresh` — no-op off main; the fresh control baseline.
- `baseline/fresh-nens500`, `exp/cr-pin-nens500` — nens=500 replicate pairs; PARKED.
- `exp/seed-rng` — determinism fixes; VERIFIED byte-identical; MERGED to main
  2026-07-02 (fast-forward to 4f96457).
- `exp/cr-pin-seeded` — seeded main + compositeR pin; paired re-test REJECTED
  (CPS +0.075, PaiCo +0.044, SCC +0.024 worse).
- `exp/vens100`, `exp/cps-neukom-target`, `exp/cps-degc` — CPS lever sweep;
  all three REJECTED (see table above).

## Project memory (copy or re-read)
Machine-local at `~/.claude/projects/-Users-nicholas-GitHub-presto-Temp12k-Composites-trial/memory/`
(`temp12k-fidelity-plan`, `-fidelity-status`, `-ground-truth-local`, `-dev-environment`).
Full plan: `~/.claude/plans/take-a-look-through-inherited-honey.md`. Everything actionable
is also summarized in this file, so a fresh session can proceed from here alone.

## nens=500 rep1 results (runs 28594449255 / 28594450895, completed 2026-07-02 ~15:10 UTC)

maxD vs published (lower better):
| method | baseline/fresh-nens500 | exp/cr-pin-nens500 | delta (pin) |
|---|---|---|---|
| SCC | 0.116 | 0.162 | +0.046 |
| DCC | 0.186 | 0.166 | -0.020 |
| GAM | 0.154 | 0.156 | +0.002 |
| CPS | 0.312 | 0.277 | -0.035 |
| PaiCo | 0.116 | 0.159 | +0.043 |

**MIXED.** Decision rule not met on rep1: PaiCo (and SCC) regressed ~+0.045; CPS/DCC
gains held direction but attenuated vs nens=100. SCC flipped sign vs nens=100 (was a
win, now a loss), so these deltas may still be noise. CPS is the only consistent
signal across nens=100 and nens=500 (pin always improves it).

## nens=500 rep2 results (runs 28600938435 / 28600940099, completed 2026-07-02 ~16:45 UTC)

maxD, rep1 / rep2 (rep1 CSVs in branch git history):
| method | baseline rep1 | baseline rep2 | pin rep1 | pin rep2 | mean delta (pin) |
|---|---|---|---|---|---|
| SCC | 0.116 | 0.139 | 0.162 | 0.115 | +0.011 |
| DCC | 0.186 | 0.144 | 0.166 | 0.167 | +0.002 |
| GAM | 0.154 | 0.155 | 0.156 | 0.155 | +0.001 |
| CPS | 0.312 | 0.375 | 0.277 | 0.376 | -0.017 |
| PaiCo | 0.116 | 0.158 | 0.159 | 0.166 | +0.026 |

### VERDICT: DO NOT PROMOTE the compositeR pin.

**Empirical nens=500 noise floor (replicate-to-replicate |diff| within a branch):**
CPS ~0.06-0.10, SCC/DCC/PaiCo ~0.04-0.05, GAM ~0.001 (deterministic).
nens=500 does NOT collapse the noise floor — run-to-run stochasticity is
structural (not tamed by ensemble size). Every pin delta, including the CPS
"win", is well inside the replicate noise. The pin is indistinguishable from
baseline. This also retroactively discredits the nens=100 single-run deltas
(CPS -0.15, DCC -0.08): with CPS wobbling ~0.1 between identical nens=500
runs, single-run comparisons cannot resolve effects of that size.

**Methodological consequence for the whole loop:** single-run A/B comparisons
on this pipeline are unreliable for all methods except GAM. Future experiments
need either (a) a pinned RNG seed in the container so comparisons are paired,
or (b) >=3 replicates per arm, comparing means. Option (a) is the cheap fix
and the recommended next infrastructure change.

`exp/cr-pin-1e3e0f2e` / `exp/cr-pin-nens500`: PARKED (not disproven, but
unresolvable at current noise). The staged Dockerfile pin stays on the branch.

## Fresh-session relaunch instructions

If the SSH session dropped and you need to restart Claude Code from scratch:
1. Open Claude Code in `/Users/nicholas/GitHub/presto-Temp12k_Composites-trial`
2. Say: "Read reproduction/ci/HANDOFF.md and continue the fidelity plan."
   Claude Code will load memories automatically and resume from this file.
3. First thing to do: run the status check above, then score if complete, or wait.
