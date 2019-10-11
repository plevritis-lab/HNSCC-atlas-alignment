#!/bin/bash
#
#SBATCH --job-name=unzip
#
#SBATCH --time=1:00:00
#SBATCH --ntasks=4
#SBATCH --wait
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
#SBATCH --requeue
#SBATCH --output=unzip.out
#SBATCH -p cmackall,owners

gunzip citeseq_36plex_$1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
gunzip citeseq_36plex_$1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
gunzip citeseq_36plex_$1/outs/filtered_feature_bc_matrix/features.tsv.gz
scp citeseq_36plex_$1/outs/filtered_feature_bc_matrix/features.tsv citeseq_36plex_$1/outs/filtered_feature_bc_matrix/genes.tsv