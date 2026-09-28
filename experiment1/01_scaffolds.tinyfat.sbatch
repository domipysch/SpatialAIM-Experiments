#!/bin/bash -l
# Stage A: one reference-scaffold build per unit (array index = line of units.tsv).
# Submit from this directory:
#   sbatch.tinyfat --array=0-251%20 01_scaffolds.tinyfat.sbatch run     (test: --array=0-2 ... test)
#SBATCH --job-name=exp1_A_scaffold
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" "${1:-run}" A
