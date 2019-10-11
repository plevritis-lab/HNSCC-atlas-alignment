#!/bin/bash
#
#SBATCH --job-name=bamtofastq
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --wait
#SBATCH --output=bamtofastq.out
#SBATCH -p cmackall

srun bamtofastq.man \
	/oak/stanford/groups/andrewg/users/zinaida/u54_data/portal.us.medgenome.com/stanford_saumyaa/$1/$1/possorted_genome_bam.bam \
	/oak/stanford/groups/andrewg/users/zinaida/u54_data/zina_fastq/$1