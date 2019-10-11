#!/bin/bash
#
#SBATCH --job-name=abseq3_7239_M
#
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=abseq3_7239_M.out
#SBATCH -p andrewg

abseq_counts \
    --bam=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/abseq_7239_M/outs/possorted_genome_bam.bam \
    --filtered_barcodes=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/cellranger_7239_M/outs/filtered_gene_bc_matrices/GRCh38 \
    --gtf=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/abseq-36plex-ref/36plex_ref.gtf \
    --rna_counts=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/cellranger_7239_M/outs/filtered_gene_bc_matrices/GRCh38 \
    --rna_info=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/cellranger_7239_M/outs/filtered_gene_bc_matrices/GRCh38 \
    --save_raw_counts \
    --output=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/combined_7239_M