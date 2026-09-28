# Shared settings for every Experiment 1 script on the NHR@FAU Tier3 clusters
# (TinyFat, frontend tinyx.nhr.fau.de). Sourced, never executed.
#
#   EXP1_BASE   everything lives here: code/, data/, envs/, and the run roots
#   EXP1_ROOT   one run root: scaffolds/, results/, logs/, summary/
#               (the one-pair test uses its own root, see README)
#   EXP1_DATA   the datasets, read in place by both stages -- the scaffold cache
#               key holds their resolved path + mtime, so never move or re-copy
#               them between stage A and stage B.

export EXP1_BASE="${EXP1_BASE:-$WORK/spatialaim_exp1}"
export EXP1_ROOT="${EXP1_ROOT:-$EXP1_BASE/run}"
export EXP1_DATA="${EXP1_DATA:-$EXP1_BASE/data/01_Datasets_Paper}"
export EXP1_CODE="${EXP1_CODE:-$EXP1_BASE/code}"
export SPATIALAIM_ENV="${SPATIALAIM_ENV:-spatialaim_env}"

# Threads per run, identical on both clusters so the timings compare.
export EXP1_THREADS="${EXP1_THREADS:-8}"

# Internet (pip, conda, GitHub) goes through the proxy on compute nodes.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  export http_proxy=http://proxy.nhr.fau.de:80
  export https_proxy=http://proxy.nhr.fau.de:80
fi

# conda from the python module; envs + package cache in $WORK, not $HOME.
# (`module` is a shell function; a non-login `bash script.sh` may lack it.)
type module >/dev/null 2>&1 || source /etc/profile
module load python >/dev/null 2>&1 || module load python
# The module sets CONDA_ENVS_PATH (read-only /apps/.../envs); conda aborts with
# MultipleKeysError when it and its alias CONDA_ENVS_DIRS are both set.
unset CONDA_ENVS_PATH
export CONDA_ENVS_DIRS="$EXP1_BASE/envs"
export CONDA_PKGS_DIRS="$EXP1_BASE/conda_pkgs"
export PIP_CACHE_DIR="$EXP1_BASE/cache/pip"
mkdir -p "$CONDA_ENVS_DIRS" "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"
_conda_base="$(conda info --base)" || { echo "config.sh: conda not usable" >&2; return 1; }
# shellcheck disable=SC1091
source "$_conda_base/etc/profile.d/conda.sh"
unset _conda_base

# Caches that would otherwise land in $HOME (the first run lost 582 tangram
# sweeps to an unwritable home). Node-local when inside a job.
_cache="${TMPDIR:-$EXP1_BASE/cache}"
export MPLCONFIGDIR="$_cache/mpl" XDG_CACHE_HOME="$_cache/xdg" \
       XDG_CONFIG_HOME="$_cache/xdg_config" NUMBA_CACHE_DIR="$_cache/numba"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$NUMBA_CACHE_DIR"
unset _cache

# Pin every threading library to the same count.
export OMP_NUM_THREADS="$EXP1_THREADS" MKL_NUM_THREADS="$EXP1_THREADS" \
       OPENBLAS_NUM_THREADS="$EXP1_THREADS" NUMBA_NUM_THREADS="$EXP1_THREADS" \
       NUMEXPR_MAX_THREADS="$EXP1_THREADS" NUMEXPR_NUM_THREADS="$EXP1_THREADS"

# Net for pandas' nullable-string default in the aligner envs (see git history).
export ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1
