# Shared settings for every Experiment 1 script on TinyFat (frontend tinyx.nhr.fau.de); sourced, never executed.
#
#   EXP1_BASE   the experiment's folder: code/, data/, envs/ and the run roots
#   EXP1_ROOT   one run root: scaffolds/, results/, logs/, summary/
#   EXP1_DATA   the datasets; never move them during a run (the scaffold cache keys on their path and mtime)

export EXP1_BASE="${EXP1_BASE:-$WORK/spatialaim_exp1}"
export EXP1_ROOT="${EXP1_ROOT:-$EXP1_BASE/run}"
export EXP1_DATA="${EXP1_DATA:-$EXP1_BASE/data/01_Datasets_Paper}"
export SPATIALAIM_ENV="${SPATIALAIM_ENV:-spatialaim_env}"

# Threads per run, the same for every job so the timings compare.
export EXP1_THREADS="${EXP1_THREADS:-8}"

# Compute nodes reach the internet only through the proxy.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  export http_proxy=http://proxy.nhr.fau.de:80
  export https_proxy=http://proxy.nhr.fau.de:80
fi

# conda from the python module (load /etc/profile first if `module` is missing).
type module >/dev/null 2>&1 || source /etc/profile
module load python >/dev/null 2>&1 || module load python
# The module's CONDA_ENVS_PATH clashes with CONDA_ENVS_DIRS below.
unset CONDA_ENVS_PATH
# Envs and package caches in $WORK, not $HOME.
export CONDA_ENVS_DIRS="$EXP1_BASE/envs"
export CONDA_PKGS_DIRS="$EXP1_BASE/conda_pkgs"
export PIP_CACHE_DIR="$EXP1_BASE/cache/pip"
mkdir -p "$CONDA_ENVS_DIRS" "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"
# Make `conda activate` available in this shell.
_conda_base="$(conda info --base)" || { echo "config.sh: conda not usable" >&2; return 1; }
# shellcheck disable=SC1091
source "$_conda_base/etc/profile.d/conda.sh"
unset _conda_base

# Tool caches off $HOME, which jobs cannot write; node-local inside a job.
_cache="${TMPDIR:-$EXP1_BASE/cache}"
export MPLCONFIGDIR="$_cache/mpl" XDG_CACHE_HOME="$_cache/xdg" \
       XDG_CONFIG_HOME="$_cache/xdg_config" NUMBA_CACHE_DIR="$_cache/numba"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$NUMBA_CACHE_DIR"
unset _cache

# Every threading library uses EXP1_THREADS.
export OMP_NUM_THREADS="$EXP1_THREADS" MKL_NUM_THREADS="$EXP1_THREADS" \
       OPENBLAS_NUM_THREADS="$EXP1_THREADS" NUMBA_NUM_THREADS="$EXP1_THREADS" \
       NUMEXPR_MAX_THREADS="$EXP1_THREADS" NUMEXPR_NUM_THREADS="$EXP1_THREADS"

# Lets anndata write pandas' nullable strings in the aligner envs.
export ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1
