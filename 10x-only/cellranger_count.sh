#!/bin/bash
#
#SBATCH --job-name=cellranger
#
#SBATCH --time=128:00:00
#SBATCH --ntasks=1
#SBATCH --wait
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --output=cellranger_count.out
#SBATCH -p cmackall

srun cellranger count \
	--id=cellranger_$1 \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/P700732_04152019/10X-only-04152019 \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-3.0.0 \
	--sample=$1 \
	--jobmode=slurm
