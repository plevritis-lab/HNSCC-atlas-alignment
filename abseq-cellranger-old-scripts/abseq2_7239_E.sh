#!/bin/bash
#
#SBATCH --job-name=abseq2_7239_E
#
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --wait
#SBATCH --output=abseq2_7239_E.out
#SBATCH -p andrewg

srun cellranger count \
	--id=abseq_7239_E \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/abseq_36plex_10xref \
	--sample=7239_E_1,7239_E_2,7239_E_3,7239_E_4 \
	--r2-length=60 \
	--jobmode=slurm