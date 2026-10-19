#!/bin/bash -l
# One K of DOT on unit 189 (pair 93, 09_human-breast-cancer, Leiden start) per array
# task, for the one run that does not fit a 24 h job; the array index is K.
# Submit from this directory, per tree (only the tree(s) without DONE):
#   sbatch.tinyfat --array=$(python 05a_run_k.py ward --levels)%20 05a_run_k.tinyfat.sbatch run ward
# then, once every K is done: python 05b_merge_k.py ward
#SBATCH --job-name=exp1_05_dot_k
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
source "$SLURM_SUBMIT_DIR/config.sh" || exit 1
case "${1:?root}" in
  run|test) export EXP1_ROOT="$EXP1_BASE/$1" ;;
  *)        export EXP1_ROOT="$1" ;;
esac
conda activate "$SPATIALAIM_ENV"
echo "=== $(date -Is)  job $SLURM_JOB_ID (array ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID})  $(hostname)"
python "$SLURM_SUBMIT_DIR/05a_run_k.py" "${2:?tree}" "$SLURM_ARRAY_TASK_ID"
