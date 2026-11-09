#!/bin/bash -l
# Demo results: one-time setup on tinyx.nhr.fau.de (run on the frontend).
#
#   SPATIALAIM_VERSION=1.5.2 bash 00_setup_server.sh   # install spatial-aim, check
#   bash 00_setup_server.sh check                       # the checks only
#
# Uses Experiment 1's envs as they are (spatialaim_env, tangram_env, tacco_env, dot_env with
# DOTr); only spatial-aim in spatialaim_env is replaced by the demo's version.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1
set -euo pipefail

SPATIALAIM_VERSION="${SPATIALAIM_VERSION:-}"
WHAT="${1:-base}"
mkdir -p "$HERE/logs"

log() { printf '\n=== %s\n' "$*"; }
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

setup_base() {
  [[ -n "$SPATIALAIM_VERSION" ]] || { echo "set SPATIALAIM_VERSION to the PyPI release of spatial-aim to use" >&2; exit 2; }
  for env in "$SPATIALAIM_ENV" tangram_env tacco_env dot_env; do
    env_exists "$env" || { echo "$env missing in $DEMO_ENVS: run Experiment 1's 00_setup_server.sh first" >&2; exit 1; }
  done
  log "Install spatial-aim==$SPATIALAIM_VERSION from PyPI"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --force-reinstall --no-deps --no-cache-dir "spatial-aim==$SPATIALAIM_VERSION"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install "spatial-aim==$SPATIALAIM_VERSION"
}

check() {
  log "spatialaim"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim --version
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python - <<'PY'
from spatialaim.reference_aligners.registry import available_reference_aligners
missing = {"tangram", "tacco", "dot"} - set(available_reference_aligners())
print("FAIL aligners missing:" if missing else "OK   all aligners available", *sorted(missing))
PY
  log "Datasets ($DEMO_DATA)"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim validate --pairs_csv "$DEMO_DATA/pairs.csv" | tail -n 5
  log "Grid"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python "$HERE/grid.py" runs
}

case "$WHAT" in
  base)  setup_base; check ;;
  check) check ;;
  *) echo "usage: $0 [base|check]" >&2; exit 2 ;;
esac
log "done"
