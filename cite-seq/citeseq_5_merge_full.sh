#!/bin/bash
#
#SBATCH --job-name=merge
#
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --wait
#SBATCH --cpus-per-task=16
#SBATCH --mem=192000
#SBATCH --requeue
#SBATCH --output=merge.out
#SBATCH -p cmackall,owners

abseq_counts \
    --bam=/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned/abseq_full_$1/outs/possorted_genome_bam.bam \
    --filtered_barcodes=/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned/citeseq_full_$1/outs/filtered_feature_bc_matrix \
    --gtf=/oak/stanford/groups/andrewg/users/zinaida/u54_data/citeseq_panels/citeseq_full_ref.gtf \
    --rna_counts=/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned/citeseq_full_$1/outs/filtered_feature_bc_matrix \
    --rna_info=/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned/citeseq_full_$1/outs/filtered_feature_bc_matrix \
    --output=/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned/combined_full_$1