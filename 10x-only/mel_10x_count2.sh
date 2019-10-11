#!/bin/bash
#
#SBATCH --job-name=cellranger
#
#SBATCH --time=128:00:00
#SBATCH --ntasks=1
#SBATCH --wait
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --output=count.out
#SBATCH -p cmackall

srun cellranger count \
	--id=cellranger_$1 \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/zina_fastq/$1/$1_MissingLibrary_1_HNJ55CCXY \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-3.0.0 \
	--jobmode=slurm
