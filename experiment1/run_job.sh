#!/bin/bash -l
# Body of every Experiment 1 array task: set up the env, then run_one.py for
# unit $SLURM_ARRAY_TASK_ID (line of units.tsv).
#
#   run_job.sh <run|test|/abs/root> A             (called by the *.sbatch scripts)
#   run_job.sh <run|test|/abs/root> B <mapper>
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1
case "${1:?root}" in
  run|test) export EXP1_ROOT="$EXP1_BASE/$1" ;;
  *)        export EXP1_ROOT="$1" ;;
esac
conda activate "$SPATIALAIM_ENV"
set -uo pipefail
mkdir -p "$EXP1_ROOT"

echo "=== $(date -Is)  job $SLURM_JOB_ID (array ${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-})  $(hostname)"
echo "root $EXP1_ROOT   stage $2 ${3:-}   cpus ${SLURM_CPUS_PER_TASK:-?}   threads $EXP1_THREADS   TMPDIR ${TMPDIR:-}"
command -v nvidia-smi >/dev/null && [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]] && nvidia-smi -L
# shellcheck disable=SC2086  # the mapper is absent for stage A
python "$HERE/run_one.py" "$2" "$SLURM_ARRAY_TASK_ID" ${3:-}
