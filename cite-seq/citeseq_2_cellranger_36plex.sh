#!/bin/bash
#
#SBATCH --job-name=cellranger
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --wait
#SBATCH --output=citeseq_count.out
#SBATCH -p cmackall

srun cellranger count \
	--id=citeseq_36plex_$1 \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-3.0.0 \
	--sample=$1_1,$1_2,$1_3,$1_4 \
	--jobmode=slurm