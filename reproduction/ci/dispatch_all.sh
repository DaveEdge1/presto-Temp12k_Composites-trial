#!/bin/bash
# Fire fidelity experiments on the fork. PREREQ: gh authed; Actions enabled on fork.
# Usage: bash dispatch_all.sh [branch ...]   (default: the three staged experiments)
set -e
REPO=${PRESTO_REPO:-nickmckay/presto-Temp12k_Composites-trial}
PAYLOAD='{"mode":"archived","compilation":"Temp12k","version":"1_0_2"}'
BR="${*:-exp/cps-t12kensemble exp/scc-median exp/cr-pin-1e3e0f2e}"
for b in $BR; do
  echo "dispatching $b"
  gh workflow run reconstruct.yml -R "$REPO" --ref "$b" \
     -f lipd_query_json="$PAYLOAD" -f unique_id="$b"
done
echo "watch: gh run list -R $REPO --workflow=reconstruct.yml"
