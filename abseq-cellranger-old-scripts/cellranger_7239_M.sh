#!/bin/bash
#
#SBATCH --job-name=cellranger_7239_M
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --wait
#SBATCH --output=cellranger_7239_M.out
#SBATCH -p andrewg

srun cellranger count \
	--id=cellranger_7239_M \
	--fastqs=/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa \
	--transcriptome=/oak/stanford/groups/andrewg/users/zinaida/refdata-cellranger-GRCh38-1.2.0 \
	--sample=7239_M_1,7239_M_2,7239_M_3,7239_M_4 \
	--jobmode=slurm