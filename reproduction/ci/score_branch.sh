#!/bin/bash
# Score a branch's CI-committed results/methods/*.csv against published (numpy-only).
# Uses `git show` (no index/working-tree side effects). Usage: score_branch.sh <branch>
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
CMP="$HERE/cmp.py"; [ -f "$CMP" ] || CMP="/Users/nicholas/.claude/jobs/586925ef/tmp/cmp.py"
b="$1"; [ -z "$b" ] && { echo "usage: score_branch.sh <branch>"; exit 1; }
cd "$ROOT"; git fetch -q origin "$b"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
got=0
for m in scc dcc gam cps paico; do
  if git show "origin/$b:results/methods/${m}_global.csv" > "$tmp/${m}.csv" 2>/dev/null; then
    got=1
    python3 "$CMP" --method "$m" --csv "$tmp/${m}.csv" 2>/dev/null \
      | grep -E "maxD|bias|spread" | tr '\n' ' ' && echo " <- $m"
  fi
done
[ "$got" = "0" ] && echo "no results/methods on origin/$b yet (CI not finished?)"
exit 0
