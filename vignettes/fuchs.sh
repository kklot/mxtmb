#!/bin/bash
#SBATCH --job-name=mx
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --time=32:00:00

export OMP_NUM_THREADS=20
srun R CMD BATCH --no-save --no-restore fit.R