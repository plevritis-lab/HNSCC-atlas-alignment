#!/bin/bash
#
#SBATCH --job-name=abseq2_7239B_M
#
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --wait
#SBATCH --output=abseq2_7239B_M.out
#SBATCH -p andrewg

srun cellranger count \
	--id=abseq_7239B_M \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/abseq_36plex_10xref \
	--sample=7239B_M_1,7239B_M_2,7239B_M_3,7239B_M_4 \
	--r2-length=60 \
	--jobmode=slurm