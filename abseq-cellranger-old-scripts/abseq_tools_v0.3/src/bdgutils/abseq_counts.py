#!/usr/bin/env python
"""
Count AbSeq reads, collapse into molecules per cell and report cell by AbSeq data tables.
"""

import sys
import os
import logging
from collections import defaultdict
import argparse
import pysam
import pandas as pd
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bdgutils.read_ref
import bdgutils.mtx_utils
import bdgutils.mtx_split
import bdgutils.demux_gene_counts
import time

LOG_FORMAT = ('[%(asctime)s] %(levelname)-8s %(name)-12s %(message)s')
logging.basicConfig(
    filename='abseq_tools.log',
    format=LOG_FORMAT,
    level=logging.DEBUG,
)
formatter = logging.Formatter(LOG_FORMAT)
handler = logging.StreamHandler()
handler.setLevel(logging.WARNING)
handler.setFormatter(formatter)
logging.getLogger().addHandler(handler)

def count_abseq_in_bam(bam, min_read, put_clabel, gtf_file, molid_tag='UB',cell_tag='CB'):
    """ count the number of reads where an AbSeq appears with a valid Cell ID and UMI
    in the bam file. 
    Return the values in a data frame where the columns are AbSeq ab-oligos and the rows are cell ids.
    If molid_tag is specified, count distinct molecule ids instead of reads.
    """
    # Store both reads and molecule counts

    def get_cell_tag(read):
        """ Get the cell_tag from the read. Return as an int if possible """
        t = read.get_tag(cell_tag)
        if t.isdigit():
            return int(t)
        else:
            return t

    bam_entries = pysam.AlignmentFile(bam, 'rb')
    contigs = [x['SN'] for x in bam_entries.header['SQ']]
    logging.info(contigs)
    for names in contigs:
        logging.info("%s", names)

    mol_counter = defaultdict(set)
    read_counter = defaultdict(int)
    readmol_counter = defaultdict(dict)
    t_counter = 0
    result_mol = defaultdict(dict)
    result_read = defaultdict(dict)
    mol_sum = defaultdict(int)
    read_sum = defaultdict(int)
    singleton_sum = defaultdict(int)
    singleton_per_cell = defaultdict(dict)

    for idx, tag in enumerate(bdgutils.read_ref.read_gtf(gtf_file)):
        # look for tag in contigs:
        gene_match = [gene for gene in contigs if gene.endswith(tag.contig)]
        if len(gene_match)==0:
            logging.info("Skipping %s", tag.contig)
            continue
        try:
            for read in bam_entries.fetch(gene_match[0], tag.start, tag.end):
                if read.is_secondary:
                    continue
                if read.get_overlap(tag.start, tag.end) < 30:
                    logging.info("length: %d, %d", tag.start, tag.end)
                    continue
                if read.has_tag(cell_tag):
                    ct = get_cell_tag(read).strip()
                else:
                    ct = 'Unknown'

                if molid_tag is not None and read.has_tag(molid_tag):
                    molid = read.get_tag(molid_tag)
                    if molid.endswith('TTTTTT'):
                        t_counter += 1
                    if molid == 0 or molid == '0' or molid == "NNNNNNNN":
                        pass
                    else:
                        # only reads with proper UMI are counted
                        read_counter[ct] += 1
                        read_sum[tag.gene_id] += 1
                        mol_counter[ct].add(molid)
                        try:
                            readmol_counter[ct][molid] += 1
                        except KeyError:
                            readmol_counter[ct][molid] = 1
                elif molid_tag is None: # count read regardless of UMI
                    read_counter[ct] += 1
                    read_sum[tag.gene_id] += 1
        except:
            continue

        for (cell, reads) in read_counter.items():
            if reads >= min_read:
                result_read[tag.gene_id][cell] = reads
        read_counter.clear()
        if molid_tag is not None:
            for (cell, mols) in mol_counter.items():
                if len(mols) > min_read:
                    result_mol[tag.gene_id][cell] = len(mols)
                    mol_sum[tag.gene_id] += len(mols)

            mol_counter.clear()

        # collect items in each tag
        for cell in readmol_counter.keys():
            singleton_per_cell[tag.gene_id][cell] = sum(i for i in readmol_counter[cell].values() if i == 1)
        readmol_counter.clear()
        singleton_sum[tag.gene_id] = sum(singleton_per_cell[tag.gene_id].values())

    with open("{}/barcodes.tsv".format(put_clabel)) as f: # need line to check it is actual cell label file later
        filtered_bcs = f.read().splitlines()
    if molid_tag is not None:
        result_mol_df = pd.DataFrame(result_mol)
        result_read_df = pd.DataFrame(result_read)
        singleton_per_cell_df = pd.DataFrame(singleton_per_cell)
        putcell_mol_df = result_mol_df.loc[result_mol_df.index.isin(filtered_bcs)]
        putcell_read_df = result_read_df.loc[result_read_df.index.isin(filtered_bcs)]
        put_singleton_per_cell_df = singleton_per_cell_df[singleton_per_cell_df.index.isin(filtered_bcs)]
    else:
        result_read_df = pd.DataFrame(result_read)
        result_mol_df = pd.DataFrame(index=result_read_df.index, columns=result_read_df.columns)
        putcell_mol_df = pd.DataFrame(index=filtered_bcs, columns=result_read_df.columns)
        singleton_per_cell_df = pd.DataFrame(singleton_per_cell)
        putcell_read_df = result_read_df.loc[result_read_df.index.isin(filtered_bcs)]
        put_singleton_per_cell_df = singleton_per_cell_df[singleton_per_cell_df.index.isin(filtered_bcs)]
    put_singleton_sum = put_singleton_per_cell_df.sum(axis=0)

    logging.info("t_counter=%s",t_counter)
    return result_read_df, result_mol_df, putcell_read_df, putcell_mol_df, singleton_sum, put_singleton_sum

def compute_stats(result_read, result_mol, putcell_read, putcell_mol, singleton_sum, put_singleton_sum, outputdir, internal):
    """Compute and output summary statistics for the AbSeq."""

    # look at putative/filtered cells
    putcell_mol_sum = []
    putcell_read_sum = []
    putcell_seq_dep = []
    num_cell_exp = []
    tot_cell_exp = 0
    ndigit = 2

    # all reads
    pAbOs = result_mol.columns.values
    mol_sum = result_mol.sum(axis=0)
    read_sum = result_read.sum(axis=0)
    tot_seq_dep = [round(x/float(y),ndigit) if y > 0 else 0 for x, y in zip(read_sum, mol_sum)]
    tot_reads = sum(read_sum)
    tot_mols = sum(mol_sum)

    # filtered cell stats
    putcell_mol_sum = putcell_mol.sum(axis=0)
    putcell_read_sum = putcell_read.sum(axis=0)
    putcell_read = putcell_read.fillna(0)
    num_cell_exp = [sum(putcell_read[pAbO].values > 0) for pAbO in pAbOs] # number of cells that has expression of each pAbO
    tot_cell_exp = sum( sum([putcell_read.sum(axis=1) > 0]) ) # number of filtered cells that express one or more pAbO
    putcell_seq_dep = [x/float(y) if y > 0 else 0 for x, y in zip(putcell_read_sum, putcell_mol_sum)]
    p_tot_reads = sum(putcell_read_sum)
    p_tot_mols = sum(putcell_mol_sum)
    num_cells = len(putcell_mol.index) # number of putative cells

    if num_cells == 0:
        mean_molpercell = np.zeros(len(pAbOs))
        tot_mean_molpercell = 0
    else:
        mean_molpercell = [x/float(num_cells) for x in putcell_mol_sum]
        tot_mean_molpercell = sum(putcell_mol_sum)/num_cells
    # output sum of filtered cells
    seq_file = os.path.join(outputdir,"stats.csv")

    filtered_bcs = putcell_read.index
    result_read = result_read.fillna(0)
    pabo_signal = result_read.loc[result_read.index.isin(filtered_bcs)]
    mean_signal = pabo_signal.mean(axis=0)
    med_signal = pabo_signal.median(axis=0)
    total_signal = result_read.loc[result_read.index.isin(filtered_bcs)].stack()
    total_mean_signal = total_signal.mean(axis=0)
    total_med_signal = total_signal.median(axis=0)

    # output statistics
    with open(seq_file, "w") as f:
        rb = csv.writer(f)
        rb.writerow(
            ['Gene', 'Num_Reads', 'Num_Molecules', 'Avg_Seq_Depth', 'Reads_Assigned_to_Cells',
             'Molecules_Assigned_to_Cells', 'Avg_Seq_Depth_Filtered_Cells',
            'Num_Filtered_Cells_with_Expression', 'Mean_Mol_Per_Cell',
             'Per_Reads_to_Filtered_Cells', 'Median_signal', 'Mean_signal'])
        for idx, tag in enumerate(pAbOs):
            list = [pAbOs[idx], read_sum[idx], mol_sum[idx], tot_seq_dep[idx],
                         putcell_read_sum[idx], putcell_mol_sum[idx], putcell_seq_dep[idx],
                         num_cell_exp[idx], mean_molpercell[idx], putcell_read_sum[idx] / float(read_sum[idx]) * 100,
                         med_signal[tag], mean_signal[tag]]
            rb.writerow([round(x,2) if i > 0 else x for i, x in enumerate(list)])
        list = ['Sum_All_Genes', tot_reads, tot_mols, round(tot_reads / float(tot_mols),ndigit),
                     p_tot_reads, p_tot_mols, p_tot_reads / float(p_tot_mols),
                     tot_cell_exp, tot_mean_molpercell, p_tot_reads / float(tot_reads) * 100,
                     total_med_signal, round(total_mean_signal,ndigit)]
        rb.writerow([round(x,2) if i > 0 else x for i, x in enumerate(list)])


def getAbSeqHistogram(df_read, df_mol, outputdir, molid_tag='UB'):
    pabo_list = df_mol.columns.values
    pabo_dir = os.path.join(outputdir, 'Histograms/')
    if not os.path.exists(pabo_dir):
        os.mkdir(pabo_dir)
    for pabo in pabo_list:
        pabo_reads = [value if value > 0 else 1 for value in df_read[pabo]]
        pabo_short = pabo.split('|')[0]
        if pabo_reads:
            plt.hist(np.log10(pabo_reads), bins=100)
            plt.xlim(xmin=0)
            plt.title(pabo_short)
            plt.ylabel('Frequency')
            plt.xlabel('log10(Reads Per Cell)')
            saveReadHist = pabo_dir + pabo_short + '_ReadsPerCell' + '.png'
            plt.savefig(saveReadHist, dpi=400)
            plt.clf()
    if molid_tag is not None:
        for pabo in pabo_list:
            pabo_mols = [value if value > 0 else 1 for value in df_mol[pabo]]
            pabo_short = pabo.split('|')[0]
            if pabo_mols:
                plt.hist(np.log10(pabo_mols), bins=100)
                plt.xlim(xmin=0)
                plt.title(pabo_short)
                plt.ylabel('Frequency')
                plt.xlabel('log10(Molecules Per Cell)')
                saveReadHist = pabo_dir + pabo_short + '_MolsPerCell' + '.png'
                plt.savefig(saveReadHist, dpi=400)
                plt.clf()

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='BAM file from a single-cell experiment with AbSeq reads where cells are identified by a CB tag.')
    parser.add_argument('--rna_counts',
                        help='Directory with matrix.mtx, barcodes.tsv, and genes.tsv files containing RNA counts from filtered cells.',required=False)
    parser.add_argument('--save_raw_counts', action='store_true')
    parser.add_argument('--rna_info', help='molecule_info.h5 file for RNA library from Cell Ranger. If provided, the script will calculate '
                                           'and save the sum of RNA and AbSeq read and molecule counts for each cell barcode in a csv file.', required=False)
    parser.add_argument('--min_read', required=False, default=0, type=int,
                        help='Minimum number of reads or molecules per protein per Cell. '
                             'Counts that do not reach this threshold will be discarded.  Used for noise reduction.')
    parser.add_argument('--output', help='Output directory.', default='./abseq_output')
    parser.add_argument('--filtered_barcodes', help='Directory with barcodes.tsv containing all the filtered cell tags')
    parser.add_argument('--gtf', help='AbSeq gtf file')
    parser.add_argument('--demux_counts_bam', help='Bam file from a single-cell experiment with BD sample tag reads', required=False)
    parser.add_argument('--demux_counts_method',
                        help="Type of Sample Tag analysis used to split cells per sample.",
                        default="noiseReduction",
                        choices=['noiseReduction',
                                 'ratioLow',
                                 'max',
                                 'cutoff_995'], required=False)
    parser.add_argument('--demux_counts_min_cell_sample_read', required=False, default=0, type=int,
                        help='Minimum number of reads per Sample Tag per Cell.  Counts that do not reach this threshold will be discarded.  Used for noise reduction.')
    parser.add_argument('--demux_counts_tag_ids', help='Comma separated list of columns in the counts file containing the sample tag data.',
                        required=False)
    args = parser.parse_args(argv)
    logging.info(("output folder = %s\n"
                  "bam file=%s\n"
                  "counts folder=%s"), args.output, args.bam, args.rna_counts)
    start_time = time.time()
    # check filtered cell barcode.tsv
    if not os.path.isfile("{}/barcodes.tsv".format(args.filtered_barcodes)):
        sys.stderr.write("Invalid directory for filtered cell barcodes.tsv")
        exit(-1)
    if not os.path.isfile(args.gtf):
        sys.stderr.write("Invalid gtf file")
    if os.path.exists(args.output):
        sys.stderr.write("{} exists\n".format(args.output))
        exit(-1)

    if args.demux_counts_bam is not None: # make sure a valid bam file is provided to count sample tags
        if not os.path.isfile(args.demux_counts_bam): # check to see if bam file exists
            sys.stderr.write("Invalid BD Sample Tag bam file (input demux_counts_bam).")
            exit(-1)
        else:
            demux_command = ['--bam', args.demux_counts_bam]
        if args.demux_counts_min_cell_sample_read is not None:
            demux_command = demux_command + ['--min_cell_sample_read', str(args.demux_counts_min_cell_sample_read)]
        if args.demux_counts_method is not None:
            demux_command = demux_command + ['--method', args.demux_counts_method]
        if args.demux_counts_tag_ids is not None:
            demux_command = demux_command + ['--tag_ids', args.demux_counts_tag_ids]

    os.mkdir(args.output)

    molid_tag='UB' # specify which tag to use for molecule count
    read_ct = False # specify whether to output read count

    # count AbSeq
    logging.info("Counting AbSeq")
    print("Begin Analysis...")
    cnt_df_read, cnt_df_mol, pcnt_df_read, pcnt_df_mol, \
    singleton_sum, put_singleton_sum  = count_abseq_in_bam(args.bam, args.min_read, args.filtered_barcodes, args.gtf, molid_tag, cell_tag='CB')
    logging.info("Count and Collapse reads.")
    print("Saving Results...")

    # save counts as csv
    if args.save_raw_counts:
        cnt_df_mol.to_csv(os.path.join(args.output, "AbSeq_MolsPerCell.csv"), index_label='barcodes')
        if read_ct:
            cnt_df_read.to_csv(os.path.join(args.output, "AbSeq_ReadsPerCell.csv"), index_label='barcodes')
    pcnt_df_mol.to_csv(os.path.join(args.output, "AbSeq_filteredMolsPerCell.csv"), index_label='barcodes')
    if read_ct:
        pcnt_df_read.to_csv(os.path.join(args.output, "AbSeq_filteredReadsPerCell.csv"), index_label='barcodes')

    if not cnt_df_read.empty:
        # output AbSeq as mtx (include all cells)
        if read_ct:
            abseq_readfile = "AbSeq_ReadsPerCell"
        abseq_molfile = "AbSeq_MolsPerCell"
        abseqrna_molfile = "AbSeqRNA_MolsPerCell"
        if args.demux_counts_bam:
            if read_ct:
                abseq_readfile = "Combined_" + abseq_readfile
            abseq_molfile = "Combined_" + abseq_molfile
            abseqrna_molfile = "Combined_" + abseqrna_molfile

        if read_ct:
            bdgutils.mtx_utils.write(pcnt_df_read.transpose(), os.path.join(args.output, abseq_readfile))
        bdgutils.mtx_utils.write(pcnt_df_mol.transpose(), os.path.join(args.output, abseq_molfile))

        if args.rna_counts is not None:
            # Get results from RNA and merge with AbSeq in mtx file
            gcm = bdgutils.mtx_utils.GeneCellMatrix.factory(args.rna_counts)  # get cell-gene matrix
            gcm.set_output_dir(args.output)  # set output directory
            # repeat for molecule counts for AbSeq
            AbSeq_mol_df = pcnt_df_mol.transpose()
            gcm.combine(AbSeq_mol_df.to_sparse(), abseqrna_molfile)

        # run demux_gene_counts if 'demux_counts' is specified'
        if args.demux_counts_bam: # separate data tables by sample tag
            if args.rna_counts is not None:
                gcm_counts = os.path.join(args.output, abseqrna_molfile)
            else: #rna gcm is not given, only output abseq
                gcm_counts = os.path.join(args.output, abseq_molfile)
            # append optional arguments
            demux_command = demux_command + ['--output', os.path.join(args.output, 'SampleTags'), '--counts', gcm_counts]
            logging.info(demux_command)
            bdgutils.demux_gene_counts.main(demux_command)

        if args.rna_info: # output total rna and protein reads /molecules for each cell label
            rna_info = bdgutils.mtx_utils.GeneCellMatrix.factory(args.rna_info)
            pabo_reads_per_cell = cnt_df_read[:-1].sum(axis=1)
            pabo_mols_per_cell = cnt_df_mol[:-1].sum(axis=1)
            protein_mol = {'protein_mols': pabo_mols_per_cell}
            protein_read = {'protein_reads': pabo_reads_per_cell}
            protein_umi_df = pd.DataFrame(data=protein_mol, index=pabo_mols_per_cell.index.values)
            protein_read_df = pd.DataFrame(data=protein_read, index=pabo_reads_per_cell.index.values)
            merged_df = pd.concat([rna_info.df, protein_read_df, protein_umi_df], axis=1)
            merged_df.fillna(0, inplace=True)
            with open("{}/barcodes.tsv".format(
                    args.filtered_barcodes)) as f:
                filtered_bcs = f.read().splitlines()
            merged_df["putative_cells"] = merged_df.index.isin(filtered_bcs) * 1
            merged_df.to_csv(args.output + '/Total_RnaAbSeqPerCell.csv', index_label='barcodes')

        print("Computing Statistics...")
        logging.info("Compute Statistics")
        # compute and output statistics
        internal = True  # if it is true, stats.csv will output more metrics
        compute_stats(cnt_df_read, cnt_df_mol, pcnt_df_read, pcnt_df_mol, singleton_sum, put_singleton_sum, args.output,
                      internal)

        print("Generating Histograms...")
        logging.info("Generating Histograms")
        # generate histograms
        getAbSeqHistogram(pcnt_df_read, pcnt_df_mol, args.output, 'UB')
        print("Done")
    else:
        logging.info("No AbSeq alignment found in bam file.")
        print("No AbSeq alignment found in bam file.")


if __name__ == "__main__":
    sys.exit(main())