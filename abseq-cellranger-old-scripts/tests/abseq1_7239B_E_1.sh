#!/bin/bash
#
#SBATCH --job-name=abseq1_7239B_E_1
#
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=abseq1_7239B_E_1.out
#SBATCH -p andrewg

cellranger count \
	--id=7239B_E_1_AbSeq \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa/7239B_E_1_Sample \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs \
	--sample=7239B_E_1 \
	--jobmode=slurm