#!/bin/bash
# Local Phase-1 runner: faithful-reproduction harness + template stack on the
# ORIGINAL v1.0.0 lpd files with real ensembles. No Docker/CI.
#
# Usage:
#   bash reproduction/localrepro/run_repro.sh cache            # build ts_all.rds (one-time, slow)
#   bash reproduction/localrepro/run_repro.sh harness <method> [nens] [seed]
#   bash reproduction/localrepro/run_repro.sh score <csv> <method>
#
# compositeR libs: rlib-1e3e0f2e (publication era) for original drivers,
# rlib-f7268c4 (container/main) for the template stack. Select via CR_LIB env:
#   CR_LIB=pub|main   (default: main)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
T12K="${T12K:-$HOME/GitHub/Temperature12k/ScientificDataAnalysis}"
LPDDIR="$T12K/lipdFilesWithEnsembles"
OUT="$HERE/cache"
mkdir -p "$OUT"

case "${CR_LIB:-main}" in
  pub)  export R_LIBS_USER="$HERE/rlib-1e3e0f2e" ;;
  main) export R_LIBS_USER="$HERE/rlib-f7268c4" ;;
  *) echo "CR_LIB must be pub|main"; exit 1 ;;
esac

cmd="${1:-}"
case "$cmd" in
  cache)
    # Build the full-TS cache + dcc slim cache (repro.R does both on first run)
    Rscript "$ROOT/reproduction/harness/repro.R" --method dcc --nens 1 \
      --lpd "$LPDDIR" --out "$OUT" --ncores 1 \
      --grid "$ROOT/reference_data/equal_area_grid_centers.csv"
    ;;
  harness)
    m="${2:?method}"; nens="${3:-500}"; seed="${4:-42}"
    case "$m" in
      dcc|scc)
        Rscript "$ROOT/reproduction/harness/repro.R" --method "$m" --nens "$nens" \
          --lpd "$LPDDIR" --out "$OUT" --ncores "${NCORES:-20}" \
          --grid "$ROOT/reference_data/equal_area_grid_centers.csv"
        ;;
      cps)
        Rscript "$ROOT/reproduction/harness/repro_cps.R" \
          --slim "$OUT/fts_cps.rds" --out "$OUT/cps_global.csv" \
          --refdata "$ROOT/reference_data" --nens "$nens" \
          --ncores "${NCORES:-20}" --seed "$seed" \
          --run-methods "$ROOT/scripts/run_methods.R"
        ;;
      paico)
        Rscript "$ROOT/reproduction/harness/repro_paico.R" \
          --slim "$OUT/fts_paico.rds" --out "$OUT/paico_global.csv" \
          --refdata "$ROOT/reference_data" --nens "$nens" \
          --ncores "${NCORES:-20}" --seed "$seed"
        ;;
      gam)
        "$HERE/venv/bin/python" "$ROOT/reproduction/harness/repro_gam.py" \
          --out "$OUT/gam_global.csv" ${GAM_ARGS:-}
        ;;
      *) echo "unknown method $m"; exit 1 ;;
    esac
    ;;
  score)
    csv="${2:?csv}"; m="${3:?method}"
    python3 "$ROOT/reproduction/ci/cmp.py" --method "$m" --csv "$csv"
    ;;
  *)
    echo "usage: run_repro.sh cache | harness <method> [nens] [seed] | score <csv> <method>"
    exit 1
    ;;
esac
