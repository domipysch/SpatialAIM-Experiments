#!/bin/bash -l
# One reference-scaffold build per pair; the array index is the PairID. Submit from this directory:
#   sbatch.tinyfat --array=$(python grid.py pairs) 01_scaffold.tinyfat.sbatch
#SBATCH --job-name=demo_scaffold
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" scaffold "$SLURM_ARRAY_TASK_ID"
