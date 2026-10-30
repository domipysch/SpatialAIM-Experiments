#!/bin/bash -l
# Experiment 1: one-time setup of all conda envs on tinyx.nhr.fau.de (run on the frontend).
#
#   SPATIALAIM_VERSION=1.5.0 bash 00_setup_server.sh   # create the envs, then run the checks
#   bash 00_setup_server.sh check                       # run the checks only
#
# Existing envs are reused; remove one first (`conda env remove -n <env>`) to rebuild it.
# spatial-aim is always reinstalled from PyPI at SPATIALAIM_VERSION.
# Tangram runs on CPU here, so tangram_env needs no GPU node.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1   # sourced before `set -u`: module/conda init read unset variables
set -euo pipefail

SPATIALAIM_VERSION="${SPATIALAIM_VERSION:-}"
WHAT="${1:-base}"
mkdir -p "$HERE/logs"   # Slurm does not create the --output folder itself

log() { printf '\n=== %s\n' "$*"; }
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# Print the installed package's reference_aligners folder (aligner env specs and wrappers).
ra_dir() {
  conda run -n "$SPATIALAIM_ENV" python -c \
    "import os, spatialaim.reference_aligners as m; print(os.path.dirname(m.__file__))"
}

# Extra packages for the runner: psutil samples memory; pandas < 3 avoids anndata's string-write error.
EXTRA_DEPS=("pandas>=2.2,<3" "psutil>=5.9")

setup_base() {
  [[ -n "$SPATIALAIM_VERSION" ]] || { echo "set SPATIALAIM_VERSION to the PyPI release of spatial-aim to pin" >&2; exit 2; }

  # spatialaim_env: Python 3.13 with spatial-aim from PyPI.
  if env_exists "$SPATIALAIM_ENV"; then
    log "$SPATIALAIM_ENV exists - reusing"
  else
    log "Create $SPATIALAIM_ENV"
    conda create -y -n "$SPATIALAIM_ENV" python=3.13 pip
  fi
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install --upgrade pip
  log "Install spatial-aim==$SPATIALAIM_VERSION from PyPI"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --force-reinstall --no-cache-dir "spatial-aim==$SPATIALAIM_VERSION" "${EXTRA_DEPS[@]}"

  # tangram_env, tacco_env, dot_env: from the env specs shipped in the package.
  local RA_DIR
  RA_DIR="$(ra_dir)"
  for name in tangram tacco dot; do
    if env_exists "${name}_env"; then
      log "${name}_env exists - reusing"
    else
      log "Create ${name}_env"
      conda env create -f "$RA_DIR/environment_${name}.yml"
    fi
  done

  # DOTr: the R package DOT's wrapper loads, installed into dot_env from GitHub.
  if dotr_loads; then
    log "DOTr R package present"
  else
    # Its R dependencies come as conda binaries: building r-maps from source fails on tinyx.
    log "Install r-fields + r-maps (conda-forge) into dot_env"
    conda install -y -n dot_env --override-channels -c conda-forge r-fields r-maps libopenblas
    # A CRAN mirror must be set, or the install only warns and still exits 0.
    log "Install the DOTr R package from GitHub (saezlab/DOT)"
    conda run --no-capture-output -n dot_env Rscript -e \
      "options(repos = c(CRAN = 'https://cloud.r-project.org')); remotes::install_github('saezlab/DOT', upgrade = 'never')"
    dotr_loads || { echo "DOTr still does not load - see the R output above" >&2; exit 1; }
    log "DOTr R package installed"
  fi
}

# True if dot_env can load DOTr and its dependency fields.
dotr_loads() {
  conda run -n dot_env Rscript -e 'suppressMessages({library(fields); library(DOTr)})' >/dev/null 2>&1
}

check() {
  # spatial-aim: installed version, and new enough for this experiment.
  log "spatialaim"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim --version
  if conda run -n "$SPATIALAIM_ENV" spatialaim run --help | grep -q no_leiden_start; then
    echo "OK   spatialaim run has --no_leiden_start"
  else
    echo "FAIL spatialaim is older than this experiment: no --no_leiden_start"
  fi
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -c "import psutil, spatialaim; print('psutil', psutil.__version__)"

  # Aligners: spatialaim finds conda and the envs, and each wrapper starts in its env.
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python - <<'PY'
from spatialaim.reference_aligners.registry import available_reference_aligners, conda_exe
print("conda:", conda_exe())
print("aligners available:", available_reference_aligners())
PY
  local RA_DIR
  RA_DIR="$(ra_dir)"
  for pair in tangram_env:run_tangram.py tacco_env:run_tacco.py dot_env:run_dot.py; do
    env="${pair%%:*}"; script="${pair##*:}"
    if env_exists "$env" && conda run -n "$env" python "$RA_DIR/$script" --help >/dev/null 2>&1; then
      echo "OK   $env runs $script"
    else
      echo "FAIL $env / $script"
    fi
  done
  if dotr_loads; then echo "OK   dot_env loads the DOTr R package"; else echo "FAIL dot_env: library(DOTr)"; fi

  # Datasets: every pair passes spatialaim's own validation.
  log "Datasets ($EXP1_DATA)"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim validate --pairs_csv "$EXP1_DATA/pairs.csv" | tail -n 5
}

case "$WHAT" in
  base)    setup_base; check ;;
  check)   check ;;
  *) echo "usage: $0 [base|check]" >&2; exit 2 ;;
esac
log "done"
