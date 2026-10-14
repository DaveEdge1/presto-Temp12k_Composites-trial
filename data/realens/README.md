# Real-ensemble bundle (build inputs, NOT committed — large)

Built by `reproduction/localrepro/build_realens_bundle.sh` from the v1.0.0 lpd
files. Contents (COPY'd into the image at build time, used when
PRESTO_REALENS=1):
- `ensemble_dcc.rds` — temp12kEnsemble+season+degC (779), 100-col ensembles
- `ensemble_cpspaico.rds` — temp12kEnsemble+season (821)
- `singlevec.json` — Temp12k+season+degC (774) single-vector proxy_ts

Regenerate when the data version changes. Gitignored (~150MB each rds).
