#!/bin/bash -l
# Body of every Experiment 1 array task: set up the env, then run_one.py with the
# remaining arguments.
#
#   run_job.sh <run|test|/abs/root> scaffold <unit>               (called by the *.sbatch scripts)
#   run_job.sh <run|test|/abs/root> sweep <unit> <mapper>
#   run_job.sh <run|test|/abs/root> k <unit> <mapper> <tree> <K>
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1
case "${1:?root}" in
  run|test) export EXP1_ROOT="$EXP1_BASE/$1" ;;
  *)        export EXP1_ROOT="$1" ;;
esac
shift
conda activate "$SPATIALAIM_ENV"
set -uo pipefail
mkdir -p "$EXP1_ROOT"

echo "=== $(date -Is)  job $SLURM_JOB_ID (array ${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-})  $(hostname)"
echo "root $EXP1_ROOT   $*   cpus ${SLURM_CPUS_PER_TASK:-?}   threads $EXP1_THREADS   TMPDIR ${TMPDIR:-}"
command -v nvidia-smi >/dev/null && [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]] && nvidia-smi -L
python "$HERE/run_one.py" "$@"
