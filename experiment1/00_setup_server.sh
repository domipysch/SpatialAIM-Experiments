#!/bin/bash -l
# Experiment 1 - one-time setup on tinyx.nhr.fau.de.
#
#   bash 00_setup_server.sh            # on the frontend: clone, spatialaim_env, tangram_env,
#                                      # tacco_env, dot_env, then the checks
#   bash 00_setup_server.sh check      # smoke checks only (frontend)
#
# Tangram runs on TinyFat's CPUs (run_tangram.py falls back to cpu without CUDA),
# so tangram_env needs no GPU node; an env built with the CUDA pytorch works too.
#
# Idempotent and non-destructive: an existing env or clone is reused; to rebuild
# one, remove it yourself first (`conda env remove -n <env>`).
#
#   SPATIALAIM_COMMIT   commit of SpatialAIM-Repository to pin (default below)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1   # before `set -u`: module/conda init read unset variables
set -euo pipefail

SPATIALAIM_COMMIT="${SPATIALAIM_COMMIT:-5a97b42375af2ae18e69a7c151f4860af2c4726b}"
REPO_URL="https://github.com/domipysch/SpatialAIM-Repository.git"
REPO="$EXP1_CODE/SpatialAIM-Repository"
RA_DIR="$REPO/src/spatialaim/reference_aligners"
WHAT="${1:-base}"
mkdir -p "$HERE/logs"   # #SBATCH --output=logs/...: Slurm does not create it

log() { printf '\n=== %s\n' "$*"; }
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# Runtime deps of spatial-aim from PyPI (pandas < 3: spatialaim's own h5ad
# writes trip anndata's nullable-string guard under pandas 3), then the pinned
# clone itself with --no-deps. psutil is for run_one.py's memory sampling.
SPATIALAIM_DEPS=(
  "numpy>=1.26" "scipy>=1.11" "pandas>=2.2,<3" "pyyaml>=6"
  "anndata>=0.12" "scanpy>=1.11" "squidpy>=1.4"
  "igraph>=0.11" "leidenalg>=0.10" "umap-learn>=0.5"
  "streamlit>=1.40" "plotly>=5.20" "kaleido>=0.2"
  "matplotlib>=3.9" "scikit-learn>=1.4" "psutil>=5.9"
)

setup_base() {
  log "Clone SpatialAIM-Repository @ $SPATIALAIM_COMMIT"
  mkdir -p "$EXP1_CODE"
  [[ -d "$REPO/.git" ]] || git clone "$REPO_URL" "$REPO"
  git -C "$REPO" fetch --quiet origin
  git -C "$REPO" checkout --quiet "$SPATIALAIM_COMMIT"
  git -C "$REPO" log -1 --oneline

  if env_exists "$SPATIALAIM_ENV"; then
    log "$SPATIALAIM_ENV exists - reusing"
  else
    log "Create $SPATIALAIM_ENV"
    conda create -y -n "$SPATIALAIM_ENV" python=3.13 pip
  fi
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install --upgrade pip
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install "${SPATIALAIM_DEPS[@]}"
  # non-editable: the installed code cannot drift if someone touches the clone
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --no-deps --force-reinstall "$REPO"

  for name in tangram tacco dot; do
    if env_exists "${name}_env"; then
      log "${name}_env exists - reusing"
    else
      log "Create ${name}_env"
      conda env create -f "$RA_DIR/environment_${name}.yml"
    fi
  done
  # The GitHub repo saezlab/DOT installs the R package DOTr (run_dot.R: library(DOTr)).
  # A non-interactive install without a CRAN mirror only warns and exits 0, so
  # set one, and verify by loading the package.
  if dotr_loads; then
    log "DOTr R package present"
  else
    # DOTr's CRAN deps fields + maps as conda binaries: building maps from source
    # on tinyx fails ("CPU ISA level is lower than required"). r-fields links
    # libopenblas.so.0 directly, which the env (libblas = mkl, via pytorch) lacks.
    log "Install r-fields + r-maps (conda-forge) into dot_env"
    conda install -y -n dot_env --override-channels -c conda-forge r-fields r-maps libopenblas
    log "Install the DOTr R package from GitHub (saezlab/DOT)"
    conda run --no-capture-output -n dot_env Rscript -e \
      "options(repos = c(CRAN = 'https://cloud.r-project.org')); remotes::install_github('saezlab/DOT', upgrade = 'never')"
    dotr_loads || { echo "DOTr still does not load - see the R output above" >&2; exit 1; }
    log "DOTr R package installed"
  fi
}

dotr_loads() {
  conda run -n dot_env Rscript -e 'suppressMessages({library(fields); library(DOTr)})' >/dev/null 2>&1
}

check() {
  log "spatialaim"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim --version
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -c "import psutil, spatialaim; print('psutil', psutil.__version__)"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python - <<'PY'
from spatialaim.reference_aligners.registry import available_reference_aligners, conda_exe
print("conda:", conda_exe())
print("aligners available:", available_reference_aligners())
PY
  for pair in tangram_env:run_tangram.py tacco_env:run_tacco.py dot_env:run_dot.py; do
    env="${pair%%:*}"; script="${pair##*:}"
    if env_exists "$env" && conda run -n "$env" python "$RA_DIR/$script" --help >/dev/null 2>&1; then
      echo "OK   $env runs $script"
    else
      echo "FAIL $env / $script"
    fi
  done
  if dotr_loads; then echo "OK   dot_env loads the DOTr R package"; else echo "FAIL dot_env: library(DOTr)"; fi
  log "Datasets ($EXP1_DATA)"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim validate --pairs_csv "$EXP1_DATA/pairs.csv" | tail -n 5
}

case "$WHAT" in
  base)    setup_base; check ;;
  check)   check ;;
  *) echo "usage: $0 [base|check]" >&2; exit 2 ;;
esac
log "done"
