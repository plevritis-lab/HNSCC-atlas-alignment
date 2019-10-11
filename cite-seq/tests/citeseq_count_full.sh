#!/bin/bash
#
#SBATCH --job-name=cellranger
#
#SBATCH --time=128:00:00
#SBATCH --ntasks=1
#SBATCH --wait
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --output=citeseq_count.out
#SBATCH -p cmackall

srun cellranger count \
	--id=citeseq_$1 \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-3.0.0 \
	--libraries=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_libraries/u54-abseq-libraries-$1.csv \
	--feature-ref=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_panels/u54-abseq-full-panel.csv \
	--jobmode=slurm
