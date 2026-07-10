# PaiCo scaling target — Neukom multi-method CFR latitude bands

These six CSVs are the 2k temperature target that PaiCo scales each zonal
composite to (per band, over the 0-1000 BP / last-millennium window; see
`scripts/paico.R`). One file per 30-deg latitude band, named by the band's
latitude range.

## Format
- `age_bp` — first column, year in BP (1950 - CE), 2000 annual rows (1949 .. -50).
- remaining columns — target **ensemble members**. `scaleComposite`/`.paico_calibrate`
  draws one member per PaiCo ensemble member and matches variance to it, which is
  what supplies PaiCo's calibration-uncertainty spread.

## Provenance
Latitude-band means of Raphael Neukom's PAGES2k-2019 climate field
reconstructions (anomalies wrt 1850-1900), i.e. the target used by the published
PaiCo pipeline (Christoph/Neukom `plotPaicoEnsemble.R`, which reads
`2k_latitude_bands_Neukom/tas_lat_bands_2k_*`). The full product is a
**600-member multi-method ensemble**: 100 members each of six reconstruction
methods (AM, CCA, CPS, DA, GraphEM, PCR).

Shipped here is a **method-balanced subsample of 198 members** (33 evenly-spaced
per method) to keep the repo footprint at ~21 MB while preserving the
between-method variance that sets PaiCo's amplitude and spread. The full
600-member target reproduces the published PaiCo to maxD 0.036; this subsample to
maxD 0.054 (both vs the single-method 100-member target's 0.112). Regenerate a
different subsample from the source files under
`~/Dropbox/Holocene GMST/2k_latitude_bands_Neukom/` if needed.

NB the earlier single-method (CPS-only) 100-member target that lived here
under-scaled the high-latitude bands and, combined with a too-wide 0-2000 BP
scaling window, produced the "too cold HTM"; both are fixed by this target + the
0-1000 window default in `paico.R`.
