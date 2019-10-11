#!/bin/bash
#
#SBATCH --job-name=mkref
#
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --output=mkref.out
#SBATCH -p cmackall

srun cellranger mkref \
	--genome=$1_reference \
	--fasta=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_panels/$1_ref.fasta \
	--genes=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_panels/$1_ref.gtf