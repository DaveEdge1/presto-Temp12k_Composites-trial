# Faithful-reproduction harness

**Reproduction-only — NOT part of the shippable Temp12k template.** This
directory holds the scripts and audit that validate this template's methods
against the published Kaufman et al. (2020) Temperature 12k results. It uses
reproduction-only conditions (the `temp12kEnsemble` record subset, the canonical
LiPD files *with* age ensembles, the published settings) that must not leak into
the user-facing template — users of the template must be free to bring their own
data outside the Temp12k compilation, with BAM age fallback.

## Layout

```
reproduction/
├── audit/                      # Fidelity scoring (anchored at 100 BP, absolute metrics)
│   ├── fidelity.py             #   per-method r, bias, RMSE, maxΔ, amp, midHol, 12ka, spread, flat
│   ├── fidelity_table.md       #   snapshot of metrics + known root causes
│   ├── overlay.py              #   6-panel overlay (mine vs published, all methods)
│   └── fidelity_overlay.png    #   the overlay figure
├── harness/                    # The reproduction pipeline (R + Python in the container)
│   ├── repro.R                 #   DCC + CPS + SCC driver (reads LiPD ensembles, composites)
│   ├── gam_dump.R              #   GAM stage 1: dump pooled cell cloud
│   ├── repro_gam.py            #   GAM stage 2: fit pygam, sample posterior, aggregate
│   ├── compare.py              #   prints harness-vs-published with the prior BAM/template line
│   ├── probe.R                 #   one-time API probe (lipdR field names, ageEnsemble dims)
│   ├── diag.R                  #   one-time NROW-mismatch / direction-field diagnostic
│   └── DCC_orig.R              #   reference: the published DCC.R driver (untouched)
└── out/                        # Committed harness outputs (for audit reproducibility)
    ├── dcc_global.csv          #   nens=50, real ageEnsemble, temp12kEnsemble subset
    ├── scc_global.csv          #   nens=50, equal-area gridding before std-over-3-5ka
    └── run_*.log               #   per-run logs (timing, record counts, filter breakdown)
```

Large caches (`ts_all.rds`, `fts_<m>.rds`, `gam_pooled.csv`) are gitignored —
they are derivable.

## Data source

The harness reads the **canonical ensemble LiPD files** from
`nickmckay/Temperature12k/ScientificDataAnalysis/lipdFilesWithEnsembles/` (698
`.lpd` files, 610 MB), which carry the real chronology ensembles the published
methods drew from. The template's runtime pipeline does NOT use these — it
downloads the lipdverse pickle (no ensembles) and BAM-simulates ages. The
divergence is intentional and the principle to preserve.

## Running it

```bash
# (1) Sparse-clone the ensemble LiPD files once
git clone --depth 1 --filter=blob:none --sparse \
    https://github.com/nickmckay/Temperature12k.git T12k
(cd T12k && git sparse-checkout set ScientificDataAnalysis/lipdFilesWithEnsembles)

# (2) Run a method via the container (presto-temp12k:local; lipdR/geoChronR/compositeR)
MSYS_NO_PATHCONV=1 docker run --rm --entrypoint bash \
  -v "$(pwd)":/repro presto-temp12k:local \
  -c "cd / && Rscript /repro/harness/repro.R --method dcc --nens 50 --out /repro/out --ncores 12"

# (3) GAM is a two-stage R→Python pipeline
docker run ... -c "Rscript /repro/harness/gam_dump.R --pool 200 --out /repro/out && \
    /opt/venv/bin/python /repro/harness/repro_gam.py --nens 50 --ncores 12"

# (4) Compare to published
python harness/compare.py --method dcc --csv out/dcc_global.csv
```

The first run builds `ts_all.rds` (~30 s to read 698 `.lpd`) and the slim
per-method cache `fts_<method>.rds` (~1 min); subsequent runs use the slim cache
directly.

## Audit method

Both my curves and the published curves are anchored to 0 at 100 yr BP (the
official NOAA convention from the readme: per ensemble member remove the full
12 ka mean, then set each method's median to 0 at 100 yr BP). Metrics are
computed on **absolute** (not recentered) anomalies, so amplitude and offset
errors are visible — unlike Pearson r / CE on recentered series, which hide
them.

See `audit/fidelity_table.md` for the latest per-method snapshot and known
root causes.

## Confirmed faithful

- **DCC** (nens=50, real `ageEnsemble`, `temp12kEnsemble` subset): spread 0.998
  (was 0.97 BAM, target 1.00), r=0.998, bias +0.035, RMSE 0.049, maxΔ 0.092.
  Real chronology ensembles fix spread perfectly. No DCC code change needed.
- **SCC** (nens=50, custom one_member_scc: bin → equal-area grid (323 cells) →
  standardize over 3-5 ka → mean cells per band): **amp 0.990** (was 1.06,
  target 1.00) — amplitude bias FIXED. midHol 0.531 vs published 0.49,
  prior 0.60. maxΔ 0.145 (prior 0.161), r=0.991, bias +0.048, RMSE 0.062.
  Residual spread 0.76 vs target 1.00 reflects an uncertainty-model difference:
  original SCC used σ=1.5 white temp noise + ±5% multiplicative age; my harness
  propagates the real value- and age-ensembles (more rigorous, narrower envelope).

## Sources of residual difference from the published (why we can't hit maxΔ→0)

1. nens=50 vs published nens=500 — sampling noise dominates maxΔ at this size
2. `compositeR` version drift — container has commit f7268c4 (~2022);
   publication-era is 1e3e0f2e (Feb 2020). `bin.R`/`spreadPaleoData`
   refactored — `minAge`, `spreadBy`, `as.numeric()` coercion differ
3. `lipdR`/`geoChronR` drift — newer extractTs returns `paleoData_values` as
   matrix; 11 records dropped due to `NROW(values)!=NROW(ageEnsemble)`
4. For SCC specifically — uncertainty-model choice (real ensembles vs σ=1.5 white)
5. Random seeds — the original didn't publish them; not recoverable
