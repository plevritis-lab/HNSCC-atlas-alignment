#!/bin/bash
#
#SBATCH --job-name=create_ref
#
#SBATCH --time=15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --output=create_ref_out
#SBATCH -p andrewg

srun create_ref --bdabseq --reflist=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/36_plex_list.txt --genome=36plex_ref