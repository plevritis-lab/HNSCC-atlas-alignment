#!/bin/bash
#
#SBATCH --job-name=slurm_test
#
#SBATCH --time=10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

srun cellranger sitecheck > sitecheck.txt
srun cellranger upload zinaida@stanford.edu sitecheck.txt
