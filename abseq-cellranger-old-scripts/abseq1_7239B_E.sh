#!/bin/bash
#
#SBATCH --job-name=abseq1_7239B_E
#
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=abseq1_7239B_E.out
#SBATCH -p andrewg

srun cellranger mkref \
	--genome=abseq_36plex_10xref \
	--fasta=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/abseq-36plex-ref/36plex_ref.fasta \
	--genes=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/abseq-36plex-ref/36plex_ref.gtf