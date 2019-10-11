#!/bin/bash
#
#SBATCH --job-name=cellranger_test
#
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --output=cellranger_test_run.out
#SBATCH -p andrewg

srun cellranger testrun --id=tiny