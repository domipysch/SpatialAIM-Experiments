# Shared settings for every demo-results script on TinyFat (frontend tinyx.nhr.fau.de); sourced, never executed.
#
#   DEMO_BASE   the demo's folder: code/, datasets/, results/, state/
#   DEMO_DATA   the demo datasets (scRNA/, ST/, pairs.csv), uploaded from Data/04_Demo_sample/datasets
#   DEMO_OUT    the finished results, laid out as the demo reads them: <linkage>/<pair_dir>/
#   DEMO_STATE  DONE/FAILED markers, run logs and the K-split sweeps; never downloaded
#   DEMO_ENVS   the conda envs; Experiment 1's, so its aligner envs are reused

export DEMO_BASE="${DEMO_BASE:-$WORK/spatialaim_demo}"
export DEMO_DATA="${DEMO_DATA:-$DEMO_BASE/datasets}"
export DEMO_OUT="${DEMO_OUT:-$DEMO_BASE/results}"
export DEMO_STATE="${DEMO_STATE:-$DEMO_BASE/state}"
export DEMO_ENVS="${DEMO_ENVS:-$WORK/spatialaim_exp1/envs}"
export SPATIALAIM_ENV="${SPATIALAIM_ENV:-spatialaim_env}"

export DEMO_THREADS="${DEMO_THREADS:-8}"

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
export CONDA_ENVS_DIRS="$DEMO_ENVS"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$WORK/spatialaim_exp1/conda_pkgs}"
export PIP_CACHE_DIR="$DEMO_BASE/cache/pip"
mkdir -p "$DEMO_BASE" "$PIP_CACHE_DIR"
_conda_base="$(conda info --base)" || { echo "config.sh: conda not usable" >&2; return 1; }
# shellcheck disable=SC1091
source "$_conda_base/etc/profile.d/conda.sh"
unset _conda_base

# Tool caches off $HOME, which jobs cannot write; node-local inside a job.
_cache="${TMPDIR:-$DEMO_BASE/cache}"
export MPLCONFIGDIR="$_cache/mpl" XDG_CACHE_HOME="$_cache/xdg" \
       XDG_CONFIG_HOME="$_cache/xdg_config" NUMBA_CACHE_DIR="$_cache/numba"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$NUMBA_CACHE_DIR"
unset _cache

export OMP_NUM_THREADS="$DEMO_THREADS" MKL_NUM_THREADS="$DEMO_THREADS" \
       OPENBLAS_NUM_THREADS="$DEMO_THREADS" NUMBA_NUM_THREADS="$DEMO_THREADS" \
       NUMEXPR_MAX_THREADS="$DEMO_THREADS" NUMEXPR_NUM_THREADS="$DEMO_THREADS"

# Lets anndata write pandas' nullable strings in the aligner envs.
export ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1
