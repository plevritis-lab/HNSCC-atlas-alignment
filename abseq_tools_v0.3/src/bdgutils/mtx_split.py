#!/usr/bin/env python
"""
Program to split a counts file into multiple files.
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as spp
import begin
import tables
import itertools
import logging
import ConfigParser


class GeneCellMatrix(object):
    @staticmethod
    def factory(mtx_file, sample=None):
        cfg = ConfigParser.ConfigParser()
        if mtx_file.endswith(".h5"):
            return HDF5_Matrix(mtx_file, "hg19", sample)
        elif os.path.isdir(mtx_file):
            return MM_Matrix(mtx_file)
        elif mtx_file.endswith("Expression_Data.st"):
            return RhapsodyST(mtx_file)
        elif mtx_file.endswith("PerCell.csv"):
            return RhapsodyDataTable(mtx_file)
        elif mtx_file.endswith(".csv"):
            return DropSeqDataTable(mtx_file)
        else:
            sys.stderr.write("Could not determine the type of file containing gene-count data.")
            return None

    def __init__(self, mtx_file, sample=None):
        self.filename = mtx_file
        self.df = self.read(sample)
        self.cells = self.df.columns
        self.genes = self.df.index
        self.sb_map = None

    def set_sample_barcode_map(self, sbm):
        """ Set the sample barcode map to be used writing out the sub-samples"""
        self.sb_map = sbm
        return self

    def set_output_dir(self, dir):
        self.output_dir = dir
        if not os.path.exists(dir):
            os.mkdir(dir)
        return self

    def read(self, sample):
        """ Read in a matrix file with gene-cell count data"""
        pass

    def write(self, sample):
        """ Write out the gene-cell count data for the cells assigned to a sample """
        pass

    def write_per_sample_data(self):
        """ Write an expression data file for each sample in the barcode_sample_map.
         The form of the output will match the form of the mtx_file.
        """
        for sample in self.sb_map.index.factorize()[1]:
            self.write(sample)

    def gene_subset(self, ids):
        " return a table (cell by tag) with all the cells, but only the genes specified by ids"
        index_ids = self.df.index.intersection(pd.Index(ids))
        return self.df.ix[index_ids].transpose().fillna(0)


class HDF5_Matrix(GeneCellMatrix):
    def __init__(self, mtx_file, genome, sample):
        self.genome = genome
        super(HDF5_Matrix, self).__init__(mtx_file, sample)

    def read(self, sample):
        with tables.open_file(self.filename, 'r') as f:
            try:
                logging.info("read h5: {}".format(f.root))
                group = f.get_node(f.root, self.genome)
                if sample:
                    group = f.get_node(group, sample)
            except:
                logging.warning("That genome ({}/{}) does not exist in this file".format(self.genome, sample))
                return None

            gene_ids = getattr(group, 'genes').read()
            gene_names = getattr(group, 'gene_names').read()
            gene_index = pd.Series(map(lambda x: "{}\t{}".format(*x), itertools.izip(gene_ids, gene_names)))
            barcodes = getattr(group, 'barcodes').read()
            data = getattr(group, 'data').read()
            indices = getattr(group, 'indices').read()
            indptr = getattr(group, 'indptr').read()
            shape = getattr(group, 'shape').read()
            matrix = spp.csc_matrix((data, indices, indptr), shape=shape)
            sdf = pd.SparseDataFrame(matrix, dtype=matrix.dtype).set_index(gene_index)
            sdf.columns = barcodes
            return sdf

    def write(self, sample):
        if sample:
            cells = self.sb_map.ix[sample, 'cell']
            if len(cells) == 0:
                return
            elif type(cells) is str:  # only 1 cell
                sdf = self.df[[cells]]
            else:
                sdf = self.df[cells]
        else:
            sdf = self.df
        with tables.open_file(self.filename, 'a') as f:
            try:
                genome_group = f.get_node(f.root, self.genome)
                sample_group = f.create_group(genome_group, sample)

            except Exception as e:
                logging.warning(repr(e))
                logging.warning("That genome ({}) does not exist in this file".format(self.genome))
                return None
            gene_lists = zip(*map(lambda x: x.split('\t'), self.genes))
            x = f.create_carray(sample_group, 'genes',
                                atom=genome_group.genes.atom,
                                shape=genome_group.genes.shape)
            x[:] = np.asarray(gene_lists[0], genome_group.genes.dtype)
            x = f.create_carray(sample_group, 'gene_names',
                                atom=genome_group.gene_names.atom,
                                shape=genome_group.gene_names.shape)
            x[:] = np.asarray(gene_lists[1], genome_group.gene_names.dtype)
            x = f.create_carray(sample_group, 'barcodes',
                                atom=genome_group.barcodes.atom,
                                shape=(len(sdf.columns),))
            x[:] = np.asarray(sdf.columns, genome_group.barcodes.dtype)
            csr = scipy.sparse.csc_matrix(sdf.to_coo(), dtype=genome_group.data.dtype)
            x = f.create_carray(sample_group, 'data', atom=genome_group.data.atom,
                                shape=csr.data.shape)
            x[:] = csr.data
            x = f.create_carray(sample_group, 'indices', atom=genome_group.indices.atom,
                                shape=csr.indices.shape)
            x[:] = csr.indices
            x = f.create_carray(sample_group, 'indptr', atom=genome_group.indptr.atom,
                                shape=csr.indptr.shape)
            x[:] = csr.indptr
            x = f.create_carray(sample_group, 'shape', atom=genome_group.shape.atom,
                                shape=(2,))
            x[:] = csr.shape
            f.close()


class MM_Matrix(GeneCellMatrix):
    def read(self, sample):
        data = scipy.io.mmread(os.path.join(self.filename, "matrix.mtx"))
        barcodes = pd.read_csv(os.path.join(self.filename, "barcodes.tsv"),
                               header=None, names=['cell'])
        genes = pd.read_csv(os.path.join(self.filename, "genes.tsv"),
                            header=None, names=['gene'])
        sdf = pd.SparseDataFrame(data).set_index(genes['gene'])
        sdf.columns = barcodes['cell']
        return sdf

    def write(self, sample):
        logging.info("output data for sample %s", sample)
        if sample:
            cells = self.sb_map.ix[sample, 'cell']
            if len(cells) == 0:
                return
            elif type(cells) is str:  # only 1 cell
                sdf = self.df[[cells]]
                sdf.columns = [cells]
            else:
                sdf = self.df[cells]
                sdf.columns = cells
        else:
            sdf = self.df
        sample_for_filename = sample.split("|")[0]
        os.mkdir(os.path.join(self.output_dir, sample_for_filename))
        pd.Series(sdf.columns).to_csv(
            os.path.join(self.output_dir, sample_for_filename, "barcodes.tsv"),
            index=False, header=False)
        pd.Series(sdf.index).to_csv(
            os.path.join(self.output_dir, sample_for_filename, "genes.tsv"),
            index=False, header=False)
        scipy.io.mmwrite(os.path.join(self.output_dir, sample_for_filename, "matrix.mtx"),
                         sdf.to_coo(),
                         comment='data for sample_for_filename: {}'.format(sample_for_filename),
                         field='integer', precision=None, symmetry=None)


class RhapsodyST(GeneCellMatrix):
    def __init__(self, mtx_file):
        self.filename = mtx_file
        self.df = self.read(None)
        self.cells = self.df.Cell_Index.unique()
        self.genes = self.df.Gene.unique()
        self.sb_map = None

    def read(self, sample):
        return pd.read_csv(self.filename, header=0, sep="\t",
                           engine='c', comment='#',
                           dtype={'DBEC_Adjusted_Molecules': np.int32,
                                  'RSEC_Reads': np.int32,
                                  'Raw_Molecules': np.int32,
                                  'RSEC_Adjusted_Molecules': np.int32,
                                  'DBEC_Reads': np.int32})

    def write(self, sample):
        sample_for_filename = sample.split("|")[0]
        if sample:
            cells = self.sb_map.ix[sample, 'cell']
            if len(cells) == 0:
                return
            sdf = self.df.set_index('Cell_Index').ix[cells].reset_index()
        else:
            sdf = self.df
        filename = os.path.join(self.output_dir,
                                "{}_{}".format(sample_for_filename, os.path.basename(self.filename)))
        sdf.to_csv(filename, sep="\t", index=False)

    def gene_subset(self, ids):
        return self.df[['Cell_Index', 'Gene', 'RSEC_Reads']].set_index(['Cell_Index', 'Gene']).unstack(1).fillna(0)[
            'RSEC_Reads'][ids]


class RhapsodyDataTable(GeneCellMatrix):
    def read(self, sample):
        return pd.read_csv(self.filename, header=0, index_col=0,
                           engine='c', comment='#').transpose()

    def write(self, sample):
        sample_for_filename = sample.split("|")[0]
        if sample:
            cells = self.sb_map.ix[sample, 'cell']
            if len(cells) == 0:
                return
            elif type(cells) is str:  # only 1 cell
                sdf = self.df[[cells]]
            else:
                sdf = self.df[cells]
        else:
            sdf = self.df
        filename = os.path.join(self.output_dir,
                                "{}_{}".format(sample_for_filename, os.path.basename(self.filename)))
        sdf.transpose().reset_index().to_csv(filename, index=False)


class DropSeqDataTable(GeneCellMatrix):
    def read(self, sample):
        df = pd.read_csv(self.filename, header=0, index_col=0,
                           engine='c', comment='#')
        if df.shape[0] > df.shape[1]:
            return df.transpose()
        else:
            return df

    def write(self, sample):
        sample_for_filename = sample.split("|")[0]
        if sample:
            cells = self.sb_map.ix[sample, 'cell']
            if len(cells) == 0:
                return
            elif type(cells) is str:  # only 1 cell
                sdf = self.df[[cells]]
            else:
                sdf = self.df[cells]
        else:
            sdf = self.df
        filename = os.path.join(self.output_dir,
                                "{}_{}".format(sample_for_filename, os.path.basename(self.filename)))
        sdf.reset_index().to_csv(filename, index=False)


@begin.start
def main(mtx_file, barcode_sample_map, output_dir='.'):
    """ Split a gene/barcode count matrix into separate files based
        on the samples in the barcode_sample_map.
    Args:
      mtx_file           : directory with matrix.mtx, barcodes.tsv, and genes.tsv files, or an h5 file, or a Rhapsody Expression_Data file, or a CSV file with Genes in columns
                           and cells in rows.
      barcode_sample_map : Pandas DataFrame or file name of a csv file containing map from cell-labels to samples.   The cell labels should match the labels in the mtx_file.
      output_dir         : if an mtx file (not h5) was specified, then this is the directory to place 
                         the output.  H5 files are modified in place; so this option would be ignored if an h5 file is
                         specified as the input.
    """
    gcm = GeneCellMatrix.factory(mtx_file)
    if isinstance(barcode_sample_map, pd.DataFrame):
        gcm.set_sample_barcode_map(barcode_sample_map.set_index('sample_tag'))
    else:
        gcm.set_sample_barcode_map(
            pd.read_csv(barcode_sample_map, header=0).set_index('sample_tag'))
    gcm.set_output_dir(output_dir).write_per_sample_data()
