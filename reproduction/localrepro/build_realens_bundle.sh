#!/bin/bash
# Build the bundled real-ensemble artifacts for the container (v1.0.0).
# Produces data/realens/{ensemble_dcc.rds, ensemble_cpspaico.rds, singlevec.json}
# from ts_all.rds (698 lpd files), column-subsampled to bound image size.
# These are BUILD INPUTS for the Docker image (COPY'd in); regenerate when the
# data version changes. Requires: reproduction/localrepro/cache/ts_all.rds and
# the slim caches (built by run_repro.sh cache + build_slim.R).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
OUT="$ROOT/data/realens"; mkdir -p "$OUT"
NCOLS="${NCOLS:-100}"
Rscript -e "
sub <- function(m,n,seed){ m<-as.matrix(m); if(ncol(m)<=n) return(m); set.seed(seed); m[,sort(sample.int(ncol(m),n)),drop=FALSE] }
subsample <- function(inp,outp){ s<-readRDS(inp)
  for(i in seq_along(s\$fTS)){ s\$fTS[[i]]\$ageEnsemble<-sub(s\$fTS[[i]]\$ageEnsemble,$NCOLS,i)
    s\$fTS[[i]]\$paleoData_values<-sub(s\$fTS[[i]]\$paleoData_values,$NCOLS,i+1e6) }
  saveRDS(s,outp,compress='xz'); cat(outp, round(file.info(outp)\$size/1e6,1),'MB\n') }
subsample('$HERE/cache/fts_dcc.rds','$OUT/ensemble_dcc.rds')
subsample('$HERE/cache/fts_cps.rds','$OUT/ensemble_cpspaico.rds')
"
cp "$HERE/cache/proxy_ts_gam.json" "$OUT/singlevec.json"
echo "bundle -> $OUT"; ls -lh "$OUT"
