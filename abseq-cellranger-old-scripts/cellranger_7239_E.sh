#!/bin/bash
#
#SBATCH --job-name=cellranger_7239_E
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --wait
#SBATCH --output=cellranger_7239_E.out
#SBATCH -p andrewg

srun cellranger count \
	--id=cellranger_7239_E \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-1.2.0 \
	--sample=7239_E_1,7239_E_2,7239_E_3,7239_E_4 \
	--jobmode=slurm