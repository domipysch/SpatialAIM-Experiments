#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Experiment 1 - server setup: create all four conda envs and install
# `spatial-aim` from TestPyPI into spatialaim_env.
#
#   bash 01_setup_envs.sh                      # full setup
#   SPATIALAIM_VERSION=0.1.4 bash 01_setup_envs.sh    # pin a TestPyPI version
#   ONLY=tangram_env bash 01_setup_envs.sh     # rebuild just one aligner env
#
# The aligner envs (tangram_env / tacco_env / dot_env) get **only their own
# aligner** - no `spatial-aim`. SpatialAIM runs their wrapper scripts by path
# (`conda run -n <env> python .../run_<aligner>.py …`), and those wrappers import
# nothing from `spatialaim`.
#
# The env specs ship inside the wheel (spatialaim/reference_aligners/environment_*.yml),
# so nothing here needs the git repository.
#
# NOTHING IS EVER DELETED: this script does not remove conda envs and does not
# touch any results directory. An env that already exists is reused as-is; to
# rebuild one, `conda env remove -n <env>` yourself first.
#
# Env vars:
#   SPATIALAIM_ENV        name of the SpatialAIM env to create           (default spatialaim_env)
#   SPATIALAIM_VERSION    TestPyPI version to pin                 (default: latest)
#   PYTHON_VERSION python for the SpatialAIM env                  (default 3.13)
#   ONLY           only build this env (spatialaim_env|tangram_env|tacco_env|dot_env)
#   SKIP_DOT_R=1   skip the DOT R package install (needs network + GitHub)
# ---------------------------------------------------------------------------
set -uo pipefail

SPATIALAIM_ENV="${SPATIALAIM_ENV:-spatialaim_env}"
PYTHON_VERSION="${PYTHON_VERSION:-3.13}"
SPATIALAIM_VERSION="${SPATIALAIM_VERSION:-}"
ONLY="${ONLY:-}"
SKIP_DOT_R="${SKIP_DOT_R:-0}"

TESTPYPI_INDEX="https://test.pypi.org/simple/"
PYPI_INDEX="https://pypi.org/simple"

# Runtime deps of spatial-aim (core + [data] extra), resolved from real PyPI.
# The package itself is then installed --no-deps from TestPyPI, so a TestPyPI
# mirror of a dependency can never be picked up by accident.
#
# pandas is capped below 3: with pandas' new string default, spatialaim's own h5ad writes
# (src/spatialaim/io.py) would hit anndata's nullable-string guard the way the aligner
# wrappers did. The wrappers are fixed; io.py is not, so keep spatialaim_env on pandas 2.
SPATIALAIM_DEPS=(
  "numpy>=1.26" "scipy>=1.11" "pandas>=2.2,<3" "pyyaml>=6"
  "anndata>=0.12" "scanpy>=1.11" "squidpy>=1.4"
  "igraph>=0.11" "leidenalg>=0.10" "umap-learn>=0.5"
  "streamlit>=1.40" "plotly>=5.20" "kaleido>=0.2"
  "matplotlib>=3.9" "scikit-learn>=1.4"
)

log()  { printf '\n=== %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

wanted() { [[ -z "$ONLY" || "$ONLY" == "$1" ]]; }

command -v conda >/dev/null 2>&1 || die "conda is not on PATH"
CONDA_BASE="$(conda info --base)" || die "conda info --base failed"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
export CONDA_EXE="${CONDA_EXE:-$CONDA_BASE/bin/conda}"

# Keep conda off $HOME: an unwritable home directory makes `conda run` abort with
# "An unexpected error has occurred" before the child process even starts. Same
# redirection as 02_run_experiment1.sh, so setup and run agree on the pkgs cache.
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$CONDA_BASE/pkgs}"
export CONDA_ENVS_DIRS="${CONDA_ENVS_DIRS:-$CONDA_BASE/envs}"
[[ -w "$CONDA_PKGS_DIRS" && -w "$CONDA_ENVS_DIRS" ]] \
  || die "conda pkgs/envs dirs are not writable ($CONDA_PKGS_DIRS, $CONDA_ENVS_DIRS)"

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# The env ymls carry their own `pip:` sections (tangram-sc, scanpy, squidpy,
# spatialdata, typst); those pip runs are launched by conda, so they cannot take a
# --no-cache-dir flag from here. The env var is how you reach them - it keeps the
# rebuild from serving anything out of a stale ~/.cache/pip and from growing it.
export PIP_NO_CACHE_DIR="${PIP_NO_CACHE_DIR:-1}"

spec="spatial-aim"
[[ -n "$SPATIALAIM_VERSION" ]] && spec="spatial-aim==$SPATIALAIM_VERSION"

# --- 1. spatialaim_env -------------------------------------------------------------
if wanted "$SPATIALAIM_ENV"; then
  if env_exists "$SPATIALAIM_ENV"; then
    log "$SPATIALAIM_ENV already exists - reusing it (remove it yourself to rebuild)"
  else
    log "Creating $SPATIALAIM_ENV (python $PYTHON_VERSION)"
    conda create -y -n "$SPATIALAIM_ENV" "python=$PYTHON_VERSION" pip || die "conda create failed"
  fi

  log "Installing dependencies from PyPI into $SPATIALAIM_ENV"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --no-cache-dir --upgrade pip || die "pip upgrade failed"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --no-cache-dir --index-url "$PYPI_INDEX" "${SPATIALAIM_DEPS[@]}" \
    || die "dependency install failed"

  # --no-cache-dir matters most here: pip would otherwise serve a wheel it already
  # cached for this version, so a re-uploaded TestPyPI build of the same version
  # would install stale code. Always fetch it fresh.
  log "Installing $spec from TestPyPI into $SPATIALAIM_ENV (--no-deps, uncached)"
  conda run --no-capture-output -n "$SPATIALAIM_ENV" python -m pip install \
    --no-cache-dir --no-deps --force-reinstall --index-url "$TESTPYPI_INDEX" "$spec" \
    || die "spatial-aim install failed"

  conda run --no-capture-output -n "$SPATIALAIM_ENV" spatialaim --version \
    || die "'spatialaim' is not runnable in $SPATIALAIM_ENV"
fi

env_exists "$SPATIALAIM_ENV" || die "$SPATIALAIM_ENV is required (it ships the aligner env specs)"

# The aligner env specs and wrapper scripts live inside the installed package.
RA_DIR="$(conda run -n "$SPATIALAIM_ENV" python -c \
  'import pathlib, spatialaim.reference_aligners as m; print(pathlib.Path(m.__file__).parent)' | tr -d '\r')"
[[ -d "$RA_DIR" ]] || die "could not locate spatialaim.reference_aligners package data"
log "Aligner env specs + wrapper scripts: $RA_DIR"

# --- 2. the three aligner envs (no spatial-aim in any of them) --------------
create_env() {   # name  yml
  local name="$1" yml="$2"
  wanted "$name" || return 0
  if env_exists "$name"; then
    log "$name already exists - reusing it (remove it yourself to rebuild)"
    return 0
  fi
  [[ -f "$yml" ]] || { warn "$yml not found - skipping $name"; return 0; }
  log "Creating $name from $yml"
  conda env create -f "$yml" || { warn "creating $name failed"; return 1; }
}

create_env tangram_env "$RA_DIR/environment_tangram.yml"
create_env tacco_env   "$RA_DIR/environment_tacco.yml"
create_env dot_env     "$RA_DIR/environment_dot.yml"

# DOT's R side: the DOT R package is not on conda, it comes from GitHub.
if wanted dot_env && env_exists dot_env && [[ "$SKIP_DOT_R" != "1" ]]; then
  if conda run -n dot_env Rscript -e 'quit(status = !requireNamespace("DOT", quietly = TRUE))' 2>/dev/null; then
    log "DOT R package already installed in dot_env"
  else
    log "Installing the DOT R package into dot_env (from GitHub)"
    conda run --no-capture-output -n dot_env \
      Rscript -e "remotes::install_github('saezlab/DOT')" \
      || warn "DOT R package install failed - the 'dot' mapper will not run"
  fi
fi

# --- 3. Verify ---------------------------------------------------------------
log "Which aligners does SpatialAIM see? (env present = usable)"
conda run --no-capture-output -n "$SPATIALAIM_ENV" python - <<'PY'
from spatialaim.reference_aligners.registry import (
    REFERENCE_ALIGNERS, available_reference_aligners, conda_exe,
)
print("conda:", conda_exe())
print("available:", available_reference_aligners())
for name, a in REFERENCE_ALIGNERS.items():
    print(f"  {name:8s} env={a.conda_env:12s} script={a.script_path}")
PY

# Each wrapper must be runnable *by path* in its own env with no spatial-aim
# installed; --help exercises its real module-level imports (torch/tacco/...).
log "Can each aligner env execute its wrapper script?"
for pair in "tangram_env:run_tangram.py" "tacco_env:run_tacco.py" "dot_env:run_dot.py"; do
  env="${pair%%:*}"; script="${pair##*:}"
  if ! env_exists "$env"; then
    warn "$env missing - '${env%%_env}' mapping will be unavailable"
    continue
  fi
  if conda run -n "$env" python "$RA_DIR/$script" --help >/dev/null 2>&1; then
    printf 'OK: %s runs %s (no spatial-aim needed)\n' "$env" "$script"
  else
    warn "$env cannot run $script - check its aligner deps"
  fi
done

log "Setup done. Next: bash 02_run_experiment1.sh (see README.md)"