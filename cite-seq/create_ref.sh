#!/bin/bash
#
#SBATCH --job-name=create_ref
#
#SBATCH --time=15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=create_ref.out
#SBATCH -p cmackall

srun create_ref \
	--bdabseq \
	--reflist=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_panels/$1_list.txt \
	--genome=$1_ref