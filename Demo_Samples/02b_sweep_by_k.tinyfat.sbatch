#!/bin/bash -l
# One K of a run too long for one 24 h job (grid.SPLIT_BY_K) per array task; the array index is K.
#   sbatch.tinyfat --array=$(python run_one.py levels <run>)%20 02b_sweep_by_k.tinyfat.sbatch <run>
# then, once every K is done: python run_one.py merge <run>
#SBATCH --job-name=demo_sweep_k
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" k "${1:?run}" "$SLURM_ARRAY_TASK_ID"
