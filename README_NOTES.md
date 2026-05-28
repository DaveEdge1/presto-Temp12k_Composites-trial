# Temperature 12k Composites

> This file is **prepended verbatim** to the auto-generated `README.md`
> by `update-readme.yml`. Everything in `README.md` below the
> `<!-- BEGIN GENERATED -->` marker is rebuilt from
> `config/user_config.yml` + `query_params.json` after every run.

**Method paper:** Kaufman, McKay, Routson, Erb & Davis (2020),
*Holocene global mean surface temperature, a multi-method reconstruction
approach*, Scientific Data 7:201, https://doi.org/10.1038/s41597-020-0530-7

## What this reconstruction does

Recreates the multi-method Holocene global mean surface temperature (GMST)
of Kaufman et al. (2020) from the Temperature 12k proxy database. Up to five
composite methods are run in one container:

| Method | Language | Approach |
|--------|----------|----------|
| **SCC** | R (compositeR) | Standard Calibrated Composite — align over 3–5 ka, equal-area gridding, native °C variance |
| **DCC** | R (compositeR) | Dynamic Calibrated Composite — iterative mean alignment, native °C variance |
| **GAM** | Python (pygam) | Per-grid-cell penalized B-spline through the age/temperature ensemble |
| **CPS** | R (compositeR) | Composite Plus Scale — z-score composite scaled to the PAGES 2k target |
| **PaiCo** | R | Pairwise Comparison — rank-based reconstruction scaled to the 2k target |

Each method composites records within six 30° latitude bands and area-weights
them into a global mean. Age uncertainty is propagated with a Banded Age Model
(symmetric 5%) generated at runtime; proxy-temperature uncertainty uses the
per-proxy-type values from the paper (Table 2). By default all five methods run
and their ensembles are **pooled into the multi-method consensus GMST** (the
paper's headline result, Fig. 3). Toggle methods off in `config/user_config.yml`
to reconstruct with a subset.

Outputs (`results/`): `reconstruction.csv` (consensus median + 5/95% + per-method
medians), `reconstruction.nc` (1D CF-NetCDF), and figures (consensus, per-method
overlay, ensemble depth).

## Pipeline

`scripts/lipd_to_ts.py` (LiPD pickle → full-fidelity `proxy_ts.json`) →
`scripts/run_methods.R` (SCC/DCC/CPS/PaiCo via compositeR) +
`scripts/gam_method.py` (GAM) → `scripts/combine_consensus.py` (pooled consensus)
→ NetCDF + figures. See [`ADAPTING.md`](ADAPTING.md) for the template scaffold.
