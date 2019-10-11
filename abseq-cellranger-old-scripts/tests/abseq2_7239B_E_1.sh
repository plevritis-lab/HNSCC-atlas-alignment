#!/bin/bash
#
#SBATCH --job-name=abseq2_7239B_E_1
#
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=abseq2_7239B_E_1.out
#SBATCH -p andrewg

abseq_counts \
    --bam=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/7239B_E_1_AbSeq/outs/possorted_genome_bam.bam \
    --filtered_barcodes=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/7239B_E_1_Sample/outs/filtered_gene_bc_matrices/GRCh38 \
    --gtf=/oak/stanford/groups/andrewg/users/zinaida/abseq-refs/36plex_ref.gtf \
    --rna_counts=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/7239B_E_1_Sample/outs/filtered_gene_bc_matrices/GRCh38 \
    --output=/oak/stanford/groups/andrewg/users/zinaida/u54_data/abseq_processed/7239B_E_1_Combined \