"""
Program to read and combine counts file from 10x cellranger and AbSeq.
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.io
import logging
import ConfigParser
import tables


class GeneCellMatrix(object):
    @staticmethod
    def factory(mtx_file):
        cfg = ConfigParser.ConfigParser()
        if os.path.isdir(mtx_file):
            return MM_Matrix(mtx_file)
        elif mtx_file.endswith("info.h5"):
            return HDF5_Info(mtx_file)
        else:
            sys.stderr.write("Could not determine the type of file containing gene-count data.")
            return None

    def __init__(self, mtx_file):
        self.filename = mtx_file
        self.df = self.read()
        self.cells = self.df.columns
        self.genes = self.df.index
        self.sb_map = None

    def set_output_dir(self, dir):
        self.output_dir = dir
        if not os.path.exists(dir):
            os.mkdir(dir)
        return self

    def read(self):
        """ Read in a matrix file with gene-cell count data"""
        pass

    def write(self):
        """ Write out the gene-cell count data """
        pass

    def combine(self, AbSeq, folder):
        " concatenate and Write out the gene-cell count data"
        pass


class MM_Matrix(GeneCellMatrix):
    def read(self):
        data = scipy.io.mmread(os.path.join(self.filename, "matrix.mtx"))
        barcodes = pd.read_csv(os.path.join(self.filename, "barcodes.tsv"),
                               header=None, names=['cell'])
        genes = pd.read_csv(os.path.join(self.filename, "genes.tsv"),
                            header=None, names=['gene'])
        sdf = pd.SparseDataFrame(data).set_index(genes['gene'])
        sdf.columns = barcodes['cell']
        return sdf

    def combine(self, AbSeq_df, folder):
        " AbSeq should be a df containing same order of cells (rows) by AbSeq barcodes (columns)"
        logging.info("combining RNA + AbSeq cell by gene data table")
        sdf = self.df
        combined_df = pd.concat([sdf, AbSeq_df])
        os.mkdir(os.path.join(self.output_dir, folder))
        pd.Series(combined_df.columns).to_csv(
            os.path.join(self.output_dir, folder, "barcodes.tsv"),
            index=False, header=False)
        # shorten gene names from 10x data table:
        genelistnew = [gene.split('\t')[1] for gene in sdf.index]
        abseqlistnew = [gene.split('|')[0]+' (AbSeq)' for gene in AbSeq_df.index]

        pd.Series(np.concatenate((genelistnew, abseqlistnew))).to_csv(
            os.path.join(self.output_dir, folder, "genes.tsv"),
            index=False, header=False)
        # pd.Series(combined_df.index).to_csv(
        #     os.path.join(self.output_dir, folder, "genes.tsv"),
        #     index=False, header=False)
        scipy.io.mmwrite(os.path.join(self.output_dir, folder, "matrix.mtx"),
                         combined_df.to_coo(),
                         comment='combined RNA+AbSeq data',
                         field='integer', precision=None, symmetry=None)


class HDF5_Info(GeneCellMatrix):

    def __init__(self, mtx_file):
        super(HDF5_Info, self).__init__(mtx_file)

    def read(self):
        with tables.open_file(self.filename, 'r') as f:
            try:
                logging.info("read h5: {}".format(f.root))
                group = f.get_node(f.root)
                print("processing")
            except:
                logging.warning("Cannot read molecule_info.h5 file")
                return None

            genome_mapped = getattr(group, 'genome').read()
            idx = np.where(genome_mapped == 0)[0]
            barcodes_num = getattr(group, 'barcode').read()[idx]
            reads = getattr(group, 'reads').read()[idx] #Number of reads that confidently mapped to this putative molecule.
            molecules = np.ones(len(idx))
            barcodes_bin = [format(bc, '032b') for bc in barcodes_num]  # convert to binary
            barcodes = [convert_to_seq(bc)+'-1' for bc in barcodes_bin]
            info = pd.DataFrame(
                {'barcode': barcodes, 'reads': reads, 'molecules': molecules})
            info.set_index('barcode', inplace=True)
            total_rna_reads = info.groupby(info.index.get_level_values(0))['reads'].apply(sum)
            total_rna_mols = info.groupby(info.index.get_level_values(0))['molecules'].apply(sum)
            combined_df = pd.concat([total_rna_reads, total_rna_mols], axis=1)
            combined_df.columns = ['rna_reads', 'rna_mols']
            combined_df.fillna(0, inplace=True)
            return combined_df


def write(abSeq_df, folder):
    logging.info("output data as mtx")
    os.mkdir(folder)
    abseqlistnew = [gene.split('|')[0] + ' (AbSeq)' for gene in abSeq_df.index]
    pd.Series(abSeq_df.columns).to_csv(
        os.path.join(folder, "barcodes.tsv"),
        index=False, header=False)
    pd.Series(abseqlistnew).to_csv(
        os.path.join(folder, "genes.tsv"),
        index=False, header=False)
    scipy.io.mmwrite(os.path.join(folder, "matrix.mtx"),
                     abSeq_df.to_sparse().to_coo(),
                     comment='data for AbSeq',
                     field='integer', precision=None, symmetry=None)


def convert_to_seq(barcodes_num):
    idx = 0
    barcodes = ''
    while idx <= len(barcodes_num)-2:
        barcodes = barcodes + bin_to_seq(barcodes_num[idx:idx+2])
        idx += 2
    return barcodes


def bin_to_seq(bin_idx):
    d = {'00':'A', '01':'C', '10':'G', '11':'T'}
    return d[bin_idx]