#!/usr/bin/env bash
# Orchestrates the Temperature 12k multi-method composite pipeline.
#
# Stages (PRESTO_STAGE env):
#   full     (default) adapter -> all methods -> consensus -> NetCDF + figures
#   methods            adapter -> the enabled method(s) -> methods/<m>_global.csv
#   combine            consensus + NetCDF + figures from existing methods/*.csv
# The CI matrix runs one `methods` job per method (PRESTO_ONLY_METHOD=<m>) in
# parallel, then a single `combine` job assembles the consensus.
#
# CI mounts:  /proxies/lipd_legacy.pkl (RO), /app/config/user_config.yml (RO), /results (RW)
# R packages (lipdR/geoChronR/compositeR) live in the renv project at / -- the R
# step runs from there so renv activates and the library resolves.

set -euo pipefail

LIPD_PICKLE="${LIPD_PICKLE:-/proxies/lipd_legacy.pkl}"
CONFIG="${PRESTO_CONFIG:-/app/config/user_config.yml}"
REFDATA="${PRESTO_REFDATA:-/app/reference_data}"
OUT="${PRESTO_OUTPUT:-/results}"
PY="${PYTHON_BIN:-/opt/venv/bin/python}"
STAGE="${PRESTO_STAGE:-full}"
ONLY="${PRESTO_ONLY_METHOD:-}"

mkdir -p "$OUT" "$OUT/methods" "$OUT/figures"
echo "[entrypoint] stage=$STAGE only_method=${ONLY:-<all>} OUT=$OUT"

# Effective config: when PRESTO_ONLY_METHOD is set (matrix job), enable just that
# one method so the container computes a single per-method global ensemble.
CONFIG_EFF="$CONFIG"
if [ -n "$ONLY" ]; then
  CONFIG_EFF="$OUT/_config_${ONLY}.yml"
  $PY -c "import yaml; c=yaml.safe_load(open('$CONFIG')) or {}; c['methods']={k:(k=='$ONLY') for k in ['scc','dcc','gam','cps','paico']}; yaml.safe_dump(c, open('$CONFIG_EFF','w'))"
fi

if [ "$STAGE" = "full" ] || [ "$STAGE" = "methods" ]; then
  echo "[entrypoint] config in use:"; cat "$CONFIG_EFF"; echo "---"

  echo "[entrypoint] Step: LiPD pickle -> proxy_ts.json"
  UNC_ARG=""
  [ -f "$REFDATA/proxy_uncertainties.yml" ] && UNC_ARG="--uncertainties $REFDATA/proxy_uncertainties.yml"
  $PY /app/scripts/lipd_to_ts.py --pickle "$LIPD_PICKLE" --out-json "$OUT/proxy_ts.json" $UNC_ARG

  echo "[entrypoint] Step: R methods (SCC/DCC/CPS/PaiCo, per config)"
  ( cd / && Rscript /app/scripts/run_methods.R \
      --ts "$OUT/proxy_ts.json" --config "$CONFIG_EFF" \
      --refdata "$REFDATA" --out-dir "$OUT/methods" )

  GAM_ON=$($PY -c "import yaml;print(bool((yaml.safe_load(open('$CONFIG_EFF')).get('methods') or {}).get('gam', True)))")
  if [ "$GAM_ON" = "True" ]; then
    echo "[entrypoint] Step: GAM (pygam)"
    $PY /app/scripts/gam_method.py --ts "$OUT/proxy_ts.json" --config "$CONFIG_EFF" \
        --grid "$REFDATA/equal_area_grid_centers.csv" --out-csv "$OUT/methods/gam_global.csv"
  fi
fi

if [ "$STAGE" = "full" ] || [ "$STAGE" = "combine" ]; then
  echo "[entrypoint] Step: consensus combine"
  $PY /app/scripts/combine_consensus.py --methods-dir "$OUT/methods" --out-csv "$OUT/reconstruction.csv"

  echo "[entrypoint] Step: CSV -> NetCDF"
  $PY /app/scripts/outputs_to_netcdf.py --in-csv "$OUT/reconstruction.csv" --out-nc "$OUT/reconstruction.nc"

  echo "[entrypoint] Step: figures"
  $PY /app/scripts/make_figures.py --in-csv "$OUT/reconstruction.csv" --out-dir "$OUT/figures"

  # results/configs.yml for the visualization / reuse UI (flatten the real config)
  $PY - "$CONFIG" "$OUT/configs.yml" <<'PY'
import sys, yaml
flat = yaml.safe_load(open(sys.argv[1])) or {}
def walk(d, p=''):
    for k, v in d.items():
        key = f'{p}.{k}' if p else k
        if isinstance(v, dict):
            yield from walk(v, key)
        else:
            yield key, v
yaml.safe_dump({'presto_config': {k: {'value': str(v)} for k, v in walk(flat)}},
               open(sys.argv[2], 'w'), default_flow_style=False)
PY
fi

echo "[entrypoint] Done (stage=$STAGE). Contents of $OUT:"
ls -lhR "$OUT" || true
