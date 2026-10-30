#!/bin/bash -l
# One reference-scaffold build per pair; the array index is the pair's Leiden unit
# (line of units.tsv). Submit from this directory:
#   sbatch.tinyfat --array=$(python grid.py array leiden)%20 01_scaffold.tinyfat.sbatch run
#   (test: --array=$(python grid.py array leiden --test) ... test)
#SBATCH --job-name=exp1_scaffold
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" "${1:-run}" scaffold "$SLURM_ARRAY_TASK_ID"
