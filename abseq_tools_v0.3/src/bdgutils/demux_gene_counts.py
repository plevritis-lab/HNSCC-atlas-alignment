#!/usr/bin/env python
"""
Count the labels on each cell and classify the cells as coming from
a particular sample.
"""

import sys
import os
import logging
from collections import defaultdict, Counter
from itertools import izip
import argparse
import pysam
import pandas as pd
import numpy as np
from scipy.stats import exponweib
from scipy.optimize import curve_fit
import bdgutils.create_ref
import bdgutils.mtx_split
import bdgutils.read_ref

LOG_FORMAT = ('[%(asctime)s] %(levelname)-8s %(name)-12s %(message)s')
logging.basicConfig(
    filename='demux_counts.log',
    format=LOG_FORMAT,
    level=logging.DEBUG,
)
formatter = logging.Formatter(LOG_FORMAT)
handler = logging.StreamHandler()
handler.setLevel(logging.WARNING)
handler.setFormatter(formatter)
logging.getLogger().addHandler(handler)


def count_tags_in_bam(bam, min_cell_sample_read, gtf_file, cell_tag='CB', molid_tag=None):
    """ count the number of reads where a sample-index-tag appears with a Cell ID
    in the bam file.
    Return the values in a data frame where the columns are sample_tags and the rows are cell ids.
    If molid_tag is specified, count distinct molecule ids instead of reads.
    """

    def get_cell_tag(read):
        """ Get the cell_tag from the read. Return as an int if possible """
        t = read.get_tag(cell_tag)
        if t.isdigit():
            return int(t)
        else:
            return t

    bam_entries = pysam.AlignmentFile(bam, 'rb')
    contigs = [x['SN'] for x in bam_entries.header['SQ']]
    logging.info("Seqname found in bam file: %s", contigs)

    # search gtf_file to determine what tags should be included
    tags_gtf = []
    if gtf_file==0: # not specified
        script_dir = os.path.dirname(__file__)
        # Try BDSampleTags:
        for tag in contigs:
            if tag.endswith('BDSampleTags'):
                tags_gtf = bdgutils.read_ref.read_gtf(os.path.join(script_dir, 'BDSampleTags.gtf'))
                logging.info("Found seqname match in default BDSampleTags-Human")
                break
            elif tag.endswith('BDSampleTagsMM'):
                tags_gtf = bdgutils.read_ref.read_gtf(os.path.join(script_dir, 'BDSampleTagsMM.gtf'))
                logging.info("Found seqname match in default BDSampleTags-Mouse")
                break
        if not tags_gtf:
            message = 'Could not find BDSampleTags or BDSampleTagsMM reference in BAM header. ' \
                      'Please specify appropriate gtf file with --tag_gtf.'
            logging.debug(contigs)
            raise Exception(message)
    else:
        tags_gtf = bdgutils.read_ref.read_gtf(gtf_file)

    mol_counter = defaultdict(set)
    read_counter = defaultdict(int)
    t_counter = 0
    result = defaultdict(dict)
    try:
        for idx, tag in enumerate(tags_gtf):
            logging.info("Reading tag %s in gtf file...", tag.contig)
            gene_match = [gene for gene in contigs if gene.endswith(tag.contig)]
            if len(gene_match) == 0:
                logging.info("Skipping %s", tag.contig)
                continue
            try:
                for read in bam_entries.fetch(gene_match[0], tag.start, tag.end):
                    if read.is_secondary:
                        continue
                    if read.get_overlap(tag.start, tag.end) < 30:
                        continue
                    if read.has_tag(cell_tag):
                        ct = get_cell_tag(read)
                    else:
                        ct = 'Unknown'
                    read_counter[ct] += 1
                    if molid_tag is not None and read.has_tag(molid_tag):
                        molid = read.get_tag(molid_tag)
                        if molid.endswith('TTTTTT'):
                            # ignore the read
                            read_counter[ct] -= 1
                            t_counter += 1
                        if molid == 0 or molid == '0' or molid == "NNNNNNNN":
                            pass
                        else:
                            mol_counter[ct].add(molid)
            except:
                continue

            if molid_tag is None:
                for (cell, reads) in read_counter.items():
                    if reads >= min_cell_sample_read:
                        result[tag.gene_id][cell] = reads
                read_counter.clear()
            else:
                for (cell, mols) in mol_counter.items():
                    if len(mols) > 10:
                        logging.debug("%s, %s", cell, mols)
                        result[tag.gene_id][cell] = len(mols)
                mol_counter.clear()
        logging.info("t_counter=%s", t_counter)
    except IndexError:
        raise Exception('Read gtf index out of range. Please check that the gtf file has valid format.')
    return pd.DataFrame(result)


def weibull_scale(data):
    """ Compute the weibull scalen parameter for a set of ratios. """
    ldata = data[data > .5]  # only look at right-hand side of distribution
    if len(ldata) > 0:
        # exponweib.fit returns a, c, loc, scale
        return exponweib.fit(ldata, floc=0)[3]
    else:
        return 0


def remove_outliers(input_list, m=4.0):
    """
    Remove outliers based on Median absolute deviation
    """
    input_list = np.array(input_list)
    deviation_from_median = np.abs(input_list - np.median(input_list))
    median_dev = np.median(deviation_from_median)
    if median_dev == 0.0:
        median_dev = 1.0
    s = deviation_from_median / median_dev
    return list(input_list[s < m])


def find_sample_tag_contig(contigs):
    """Return the first contig found that references the sample tag chromosome name"""
    for contig in contigs:
        if contig.endswith('BDSampleTags'):
            return contig
    else:
        message = 'Could not find BDSampleTags reference in BAM header.'
        logging.debug(contigs)
        raise Exception(message)


class SampleChooser(object):
    @staticmethod
    def factory(name, data):
        if name == 'cutoff_995':
            return CutoffChooser(data)
        if name == 'ratioLow':
            return RatioChooser(data, .35)
        if name == 'ratioHigh':
            return RatioChooser(data, .75)
        if name == 'max':
            return MaxChooser(data, 10)
        if name == 'noiseReduction':
            return NoiseReductionChooser(data, .75, 30)

        return None

    def choose_label(self, label_vector):
        result = []
        for i in range(len(label_vector)):
            if label_vector[i]:
                result.append(self.tags[i])  # collect name of column
        if len(result) == 1:
            return result[0]
        elif len(result) == 0:
            return "Undetermined"
        else:
            return "Multiplet"

    def assign_cells():
        pass


class CutoffChooser(SampleChooser):
    """
        This class sets a cutoff based on fitting a normal curve to the log-transformed tag counts,
        and selecting a cutoff at 2.33 std deviations from the mean.

        In practice this tends to choose very low values for the cutoff and overcalls multiplets.

    """

    def __init__(self, data):
        # drop any sample tags that don't have at least 1000 reads.
        self.df = data.fillna(0).drop([col for col in data.columns if data[col].sum() <= 1000], axis=1)
        self.levels = []
        self.assignments = []
        self.name = "Gaussian Cutoff Method"
        logging.info("%s, columns = %s", self.name, self.df.columns)
        self.tags = self.df._get_numeric_data().columns

    @staticmethod
    def gauss(x, *p):
        A, mu, sigma = p
        return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

    @staticmethod
    def get_cutoff(data, p0):
        """
        fit a gaussian curve on the log-transformed data and find a cutoff value
        """

        ldata = data.map(np.log1p)
        if sum(ldata) == 0.0:
            return None, None

        # if curve_fit fails, try skipping more bins to get rid of
        # noise at the low-end of the distribution.
        for skip in range(1, 10):
            hist, bin_edges = np.histogram(ldata.values, bins=20)
            bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
            # Define model function to be used to fit to the data above:
            # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
            if p0 is None:
                x = np.arange(len(ldata))
                mean = sum(x * ldata) / sum(ldata)
                sigma = np.sqrt(sum(ldata * (x - mean) ** 2) / sum(ldata))
                p0 = [max(ldata), mean, sigma]
            try:
                (A, mean, sigma), var_matrix = curve_fit(CutoffChooser.gauss,
                                                         bin_centres[skip:],
                                                         hist[skip:], p0=p0)
                logging.info("skip, mean, sigma = %d, %6.4f, %6.4f",
                             skip, mean, sigma)
                if mean - (2.33 * abs(sigma)) >= .7:  # keep going if this doesn't remove the singletons.
                    return np.exp(mean - (2.33 * abs(sigma))) - 1, [A, mean, sigma]
            except Exception as e:
                logging.warning("failed with skip= %s : %s", skip, e)
        return None, None

    def assign_cells(self):
        """ add a cluster annotatation to each cell.
        Returns:
          updates the object with an variable 'assignments' containing the list of samples corresponding to each cell.
          The assignments are either Undetermined, Multiplet, or the column name for the single column that is above
          a computed cutoff value for the column.
        """
        # find thresholds for each column
        # levels is a list floats indicating the highest count still considered noise for each column

        mx = self.df.max(axis=1)
        p0 = CutoffChooser.get_cutoff(mx, None)[1]
        self.levels = [(col,
                        CutoffChooser.get_cutoff(self.df[col], p0)[0]) \
                       for col in self.tags]
        logging.info("%s", self.levels)
        labeled = [self.df[col] > x for (col, x) in self.levels]
        self.assignments = [self.choose_label(x) for x in izip(*labeled)]
        logging.info("sample_tag count = %s", Counter(self.assignments))
        assert len(self.assignments) > 0
        return self.assignments


class NoiseReductionChooser(SampleChooser):
    """
        This class sets a cutoff based on finding a separation between noise and signal values for each sample.  It attempts to remove
        the typical noise values from each sample before computing the assignment.   If there is insufficient separation it uses the highest
        signal relative to the noise value.
    """

    def __init__(self, data, ratio, min_count_for_singleton):
        # drop any sample tags that don't have at least 1000 reads.
        self.df = data.fillna(0).drop([col for col in data.columns if data[col].sum() <= 1000], axis=1)
        self.ratio = ratio  # good ratio threshhold for computing noise.
        self.minimum_count_for_singleton = min_count_for_singleton
        self.assignments = []
        self.levels = None
        self.status = {}
        self.name = "NoiseReduction Method: {:4.3f}".format(ratio)
        self.tags = self.df._get_numeric_data().columns
        logging.info("%s, columns = %s", self.name, self.df.columns)

    def _characterize_noise(self):
        " compute trend_noise by looking at the relationship between good counts associated with a cell and the noise for each cell "
        total_count_x = self.df[self.tags].sum(axis=1)
        max_vals = self.df[self.tags].max(axis=1)
        noise_count_y = total_count_x - max_vals
        self.total_count_x = total_count_x
        self.good_cells = ((max_vals / total_count_x.astype(float)) > self.ratio) & (
        max_vals >= self.minimum_count_for_singleton)
        if self.good_cells.any():
            logging.info("good cells = %d, %s", len(self.good_cells), self.good_cells)
            trend_z = np.polyfit(total_count_x[self.good_cells], noise_count_y[self.good_cells], 1)
            logging.debug("%s\n\n%s", total_count_x[self.good_cells], noise_count_y[self.good_cells])
            self.trend_noise = np.poly1d(trend_z)
            return True
        else:
            self.levels = {x: 1 for x in self.tags}
            self.status = {x: 'Warning - No unambiguous cells detected' for x in self.tags}
            return False

    def _determine_cutoff_per_tag(self):
        """
        Determine the minimum read-count required make a call for each sample.   The try to find the lowest value that's clearly separated from
        the noise.
        """
        self.max_tags = self.df[self.tags].apply(pd.Series.idxmax, axis=1)
        self.tag_to_good_counts = {}
        tag_to_noise_counts = {}
        for tag in self.tags:
            self.tag_to_good_counts[tag] = self.df[tag][(self.max_tags == tag) & self.good_cells]
            tag_to_noise_counts[tag] = self.df[tag][(self.max_tags != tag) & (self.df[tag] > 0) & self.good_cells]
            if len(tag_to_noise_counts[tag]) > 0:
                tag_to_noise_counts[tag] = remove_outliers(tag_to_noise_counts[tag], m=12.0)

        tag_to_min_good = {}

        for tag, cnts in self.tag_to_good_counts.items():
            if len(cnts) > 0:
                # Remove outliers from each tag_to_good_counts[tag]
                # Only consider the lowest 100 counts
                good_counts_low_segment = sorted(cnts)[:100]
                logging.info("gcls: %s", good_counts_low_segment)
                clean_good_counts_low_segment = remove_outliers(good_counts_low_segment, m=4.0)
                logging.info("cgcls: %s", clean_good_counts_low_segment)
                # Allow some flexibility for the minimum good count to be below the determined min good count
                #   e.g. 80% with this:  (np.min(tag_to_good_counts[tag]) * 0.8)
                # Lowering this value will increase the number of multiplets
                tag_to_min_good[tag] = np.min(clean_good_counts_low_segment) * 0.75
                logging.info("tag_to_min_good[%s]:%s", tag, tag_to_min_good[tag])
                logging.info("max(tag_to_noise_counts[%s]:%s", tag, max(tag_to_noise_counts[tag]))
                #  ... But make sure we have separation of good count and noise.  So min_good must higher than outlier cleaned noise counts
                if len(tag_to_noise_counts[tag]) > 0:
                    if tag_to_min_good[tag] < max(tag_to_noise_counts[tag]):
                        logging.info('Initial min_good (%s) less than noise', tag_to_min_good[tag])
                        # make the min_good 50% higher than the maximum cleaned noise count
                        logging.info("tag_to_noise_counts[%s]:%s", tag, tag_to_noise_counts[tag])
                        tag_to_min_good[tag] = max(tag_to_noise_counts[tag]) * 1.5
                        self.status[tag] = 'Warning - Signal and Noise too close'
                        logging.info('Warning - Signal and Noise too close for %s', tag)
                    else:
                        self.status[tag] = 'Good'
                        logging.info('Sample Tag Signal for %s -- good', tag)
            else:
                # There were no good singleton counts for this tag.
                tag_to_min_good[tag] = None

        self.levels = tag_to_min_good
        for tag in tag_to_min_good:
            tag_to_noise_counts[tag] = [noise_count for noise_count in tag_to_noise_counts[tag] if
                                        noise_count < tag_to_min_good[tag]]

        # Count up the noise counts for each tag so that we know their contribution to the noise relative to the overall noise
        tag_to_noise_total = {tag: sum(tag_to_noise_counts[tag]) for tag in tag_to_noise_counts}
        total_noise = sum(tag_to_noise_total.values())
        if total_noise > 0:
            self.tag_to_percent_noise = {tag: tag_to_noise_total[tag] / float(total_noise) for tag in
                                         tag_to_noise_total}
        else:
            self.tag_to_percent_noise = {tag: 0 for tag in tag_to_noise_total}

    def _remove_noise(self):
        """
        Remove expected noise from each measurement by fitting the trend_noise curve, and using the percent of noise from each tag
        to scale that value appropriately for each sample.

        This method modifies the objects df component.

        """
        expected_noise = self.trend_noise(self.total_count_x)
        noise_per_tag = pd.DataFrame({tag: self.tag_to_percent_noise[tag] * expected_noise for tag in self.tags},
                                     index=self.df.index).round(0).astype(int)
        logging.info("noise_per_tag: %s", noise_per_tag)
        for t in self.tags:
            self.df[t] = np.maximum(self.df[t] - noise_per_tag[t], 0)

    def _choose_best_assignment_for_unassigned_cells(self):
        """
        For cells that don't have a value tag over the minimum good cutoff we choose the highest value relative to the noise values.

        This method updates the assignments component.
        """
        min_normed = .2
        min_count = 5
        normed = self.df[self.tags].apply(lambda x: (x[x >= min_count] / self.levels[x.name]), axis=0)
        logging.debug("normed: %s", normed)
        normed_max = normed.max(axis=1)
        logging.debug("normed_max: %s", normed_max)
        labeled = [(normed[col] - normed_max == 0) & (normed[col] > min_normed) for col in self.tags]
        logging.debug("choosing best assignment for unassigned: %s", labeled)
        assignment = pd.Series([self.choose_label(x) for x in izip(*labeled)], index=normed.index)
        logging.debug("assignment: %s", assignment)
        unassigned = (self.assignments == 'Undetermined')
        logging.info("unassigned = %s", unassigned)
        self.assignments[unassigned] = assignment[unassigned]
        self.assignments.fillna('Undetermined', inplace=True)

    def assign_cells(self):
        """ add a cluster annotation to each cell.
        Returns:
          updates the object with an variable 'assignments' containing the list of samples corresponding to each cell.
          The assignments are either Undetermined, Multiplet, or the column name for the single column that is above
          a computed cutoff value for the column.
        """
        if self._characterize_noise():
            self._determine_cutoff_per_tag()
            self.tags = [t for t in self.tags if self.levels[t] is not None]
            self._remove_noise()
            labeled = [self.df[t] >= self.levels[t] for t in self.tags]
            logging.debug("%s", labeled)
            self.assignments = pd.Series([self.choose_label(x) for x in izip(*labeled)], index=self.df.index)
        else:
            self.assignments = pd.Series(['Undetermined' for x in self.df.index], index=self.df.index)

        self._choose_best_assignment_for_unassigned_cells()
        logging.info("sample_tag count = %s", Counter(self.assignments))
        assert len(self.assignments) > 0
        return self.assignments


class RatioChooser(SampleChooser):
    """
    This class chooses the sample based on the ratio of the (maximum tag read count)/(sum of all tag read counts tags).  All the samples above a minimum ratio are considered valid samples for the tag, and cells with multiple assignments are called Multiplet.

    Before computing the ratio, the data is normalized so that the total # reads per sample is equal.

    Only samples with a 'passing' status are considered. The status is determined by computing fitting a Weibull curve, and testing
    that the scale parameter is > .5.

    The method is highly sensitive to the minimum ratio used.  The ratioLow method uses .35 which tends to overcall Multiplets.
    """

    def __init__(self, data, ratio):
        # drop any sample tags that don't have at least 1000 reads.
        self.df = data.fillna(0).drop([col for col in data.columns if data[col].sum() <= 1000], axis=1)
        self.ratio = ratio
        self.assignments = []
        self.levels = None
        self.name = "Ratio Method: {:4.3f}".format(ratio)
        logging.info("%s, columns = %s", self.name, self.df.columns)
        self.tags = self.df._get_numeric_data().columns

    def assign_status(self):
        sums = self.df.sum(axis=1)
        norm_factor = self.df.sum(axis=0).sum() / self.df.sum(axis=0)
        self.weib = pd.Series([weibull_scale(self.df[t] * norm_factor[t] / sums) for t in self.tags],
                              index=self.tags)
        logging.info("weibull values: %s", self.weib)
        self.status = self.weib > .5
        logging.info("weibull status: %s", self.status)

    def tags_that_fail(self):
        return [x[0] for x in izip(self.tags, self.status) if x[1] == False]

    def assign_cells(self):
        """ add a cluster annotation to each cell.
        Returns:
          updates the object with an variable 'assignments' containing the list of samples corresponding to each cell.
          The assignments are either Undetermined, Multiplet, or the column name for the single column that is above
          a computed cutoff value for the column.
        """
        self.assign_status()
        self.df = self.df.drop(self.tags_that_fail(), axis=1)
        self.tags = self.df._get_numeric_data().columns
        sums = self.df.sum(axis=1)
        norm_factor = sum(self.df.sum(axis=0)) / self.df.sum(axis=0)
        logging.info("norm_factors = %s", norm_factor)
        labeled = [(self.df[col] * norm_factor[col]) / sums > self.ratio for col in self.tags]
        self.assignments = [self.choose_label(x) for x in izip(*labeled)]
        logging.info("sample_tag count = %s", Counter(self.assignments))
        assert len(self.assignments) > 0
        return self.assignments


class MaxChooser(SampleChooser):
    """
    Classify cells based on the maximum read count.
    The data is not normalized or filtered (except to remove samples with very low total counts.)

    Almost no multiplets are called.   A multiplet will only be called when a cell has the same highest tag count from two samples.
    """

    def __init__(self, data, min_callable):
        self.df = data.fillna(0).drop([col for col in data.columns if data[col].sum() <= 1000], axis=1)
        self.min_callable = min_callable
        self.assignments = []
        self.levels = None
        self.name = "Max Value Method"
        logging.info("%s, columns = %s", self.name, self.df.columns)
        self.tags = self.df._get_numeric_data().columns

    def assign_cells(self):
        """ add a cluster annotatation to each cell.
        Returns:
          updates the object with an variable 'assignments' containing the list of samples corresponding to each cell.
          The assignments are either Undetermined, Multiplet, or the column name for the single column that is above
          a computed cutoff value for the column.
        """
        maxes = self.df.max(axis=1)
        labeled = [self.df[col] - maxes == 0 for col in self.tags]
        self.assignments = [self.choose_label(x) for x in izip(*labeled)]
        logging.info("sample_tag count = %s", Counter(self.assignments))
        assert len(self.assignments) > 0
        return self.assignments


def compute_stats(all_df, filtered_cntdf, chooser):
    """Compute summary statistics for the sample tags."""

    def getV(gc, v, c):
        if v in gc.index and c in gc.columns:
            return gc.ix[v][c]
        else:
            return 0

    def zero_div(x,y):
        return x / y if y else 0

    method = chooser.keys()[0]
    experiment_stats = dict()
    experiment_stats['num_samples'] = [len(set(chooser[method].assignments).difference(set(['Undetermined', 'Multiplet'])))]
    experiment_stats['num_cells'] = [len(all_df)]
    experiment_stats['num_filtered_cells'] = [len(chooser[method].df)]
    filtered_cntdf.fillna(0, inplace=True)
    groups = [filtered_cntdf.reset_index().groupby([method])]

    samples = set()
    for g in groups:
        samples.update(g.count().index)
        sample_tags = list(g.count().columns)[1:]
        tag_calls = list(g.groups.keys())
        g = g

    for s in samples:
        experiment_stats[s] = [getV(g.count(), s, 'cell') for g in groups]

    if 'Undetermined' not in experiment_stats:
        experiment_stats['Undetermined'] = [0 for g in groups]
    if 'Multiplet' not in experiment_stats:
        experiment_stats['Multiplet'] = [0 for g in groups]

    if float(experiment_stats['Undetermined'][0])/experiment_stats['num_filtered_cells'][0] > 0.1:
        message = "Warning: High percentage ({}%) of cells have no sample tag assignment".format(
            round(float(experiment_stats['Undetermined'][0])/experiment_stats['num_filtered_cells'][0]*100, 2))
        logging.info(message)
        print(message)

    experiment_stats['num_uniquely_assigned'] = [x[0] - x[1] - x[2] for x in
                                                 izip(experiment_stats['num_filtered_cells'],
                                                      experiment_stats['Multiplet'],
                                                      experiment_stats['Undetermined'])]

    SNR = defaultdict(dict)
    for s in sample_tags:
        if s in tag_calls:
            SNR['mean_signal'][s] = g[s].mean()[s]
            SNR['median_signal'][s] = g[s].median()[s]
            SNR['mean_noise'][s] = filtered_cntdf[filtered_cntdf[method].isin([st for st in sample_tags if st != s])][s].mean()
            SNR['median_noise'][s] = filtered_cntdf[filtered_cntdf[method].isin([st for st in sample_tags if st != s])][s].median()
            SNR['mean_SNR'][s] = zero_div(SNR['mean_signal'][s], SNR['mean_noise'][s])
            SNR['median_SNR'][s] = zero_div(SNR['median_signal'][s], SNR['median_noise'][s])
            if SNR['mean_SNR'][s] < 10.0:
                logging.info("Low signal to noise ratio in %s = %0.2f", s, round(SNR['mean_SNR'][s],2))
                print("Warning: Low signal to noise ratio in {} = {}".format(s, round(SNR['mean_SNR'][s],2)))

    exp_stats = pd.DataFrame(experiment_stats, index=[method])

    def count_non_singletons(x):
        return x[x > 1].count()

    def mean_non_singletons(x):
        return np.mean(x[x > 1])

    sample_tag_stats = chooser[method].df.agg([np.sum, np.mean],
                                              axis=0).transpose()
    sample_tag_stats.rename(columns={'sum':'Total_Reads_in_Filtered_Cells', 'mean':'Mean_Reads_per_Filtered_Cell'}, inplace=True )
    sample_tag_stats.index.rename('sample_tag', inplace=True)
    SNR_dt = pd.DataFrame(SNR)
    SNR_dt.index.rename('sample_tag', inplace=True)
    sample_tag_stats = pd.merge(sample_tag_stats, SNR_dt, how='left', left_index=True,
                                right_index=True)

    all_sample_tag_stats = all_df.agg([np.sum, np.mean],
                                      axis=0).transpose()
    all_sample_tag_stats.rename(columns={'sum':'Total_Reads', 'mean':'Mean_Reads_per_Barcode'}, inplace=True)
    sample_tag_stats = pd.merge(all_sample_tag_stats, sample_tag_stats, how='left', left_index=True, right_index=True)

    logging.debug("sample_tag_stats.columns= %s", sample_tag_stats.columns)
    logging.debug("sample tag stats = %s", sample_tag_stats)
    logging.debug("%s = %s", method, chooser[method].levels)

    return (exp_stats, sample_tag_stats)


def write_summary_stats(output, chooser, exp_stats, sample_tag_stats):
    """ Output the computed stats. """

    output.write("## Summary Stats for Sample Demultiplexing ##\n")
    exp_stats.to_csv(output)

    if sample_tag_stats is not None:
        output.write("\n## Sample Tag Stats (read counts)##\n")
        sample_tag_stats.round(2).to_csv(output)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='BAM file from a single-cell experiment where cells are identified by a CB tag.',
                        required=False)
    parser.add_argument('--tag_counts', help='csv file containing a data table with sample tag counts for each cell label.', required=False)
    parser.add_argument('--tag_ids',
                        help='Comma separated list of columns in the counts file containing the sample tag data.',
                        required=False)
    parser.add_argument('--barcodes', help='List of cell barcodes produced by mRNA data with no header. '
                                           'Only needed if mRNA data table is not provided in --counts', required=False)
    parser.add_argument('--counts',
                        help='Directory with matrix.mtx, barcodes.tsv, and genes.tsv files, or an h5 file, or a '
                             'Rhapsody Expression_Data file, or a CSV file with Genes in columns and cells in rows.',
                        required=False)
    parser.add_argument('--method',
                        help="Type of Sample Tag analysis used to split cells per sample.",
                        default="noiseReduction",
                        choices=['noiseReduction',
                                 'ratioLow',
                                 'max',
                                 'cutoff_995'])
    parser.add_argument('--tag_gtf', help='Specify gtf file used, if it differs from the default human or mouse gtf '
                                          'generated by multiplexing_tools', default=0)
    parser.add_argument('--save_raw_counts', action='store_true')
    parser.add_argument('--min_cell_sample_read', required=False, default=0, type=int,
                        help='Minimum number of reads per Sample Tag per Cell.  '
                             'Counts that do not reach this threshold will be discarded.  Used for noise reduction.')
    parser.add_argument('--output', help='Output directory.', default='./demux_output')
    args = parser.parse_args(argv)
    logging.info(("output folder = %s\n"
                  "bam file=%s\n"
                  "counts folder=%s"), args.output, args.bam, args.counts)

    # check that tag_gtf is valid
    if not args.tag_gtf==0:
        if (not os.path.isfile(args.tag_gtf)) or (not args.tag_gtf.endswith('.gtf')):
            sys.stderr.write("Invalid gtf file (--tag_gtf)\n")
            exit(-1)

    # check that either a list of cell barcodes or counts
    if args.counts is not None:
        if not os.path.exists(args.counts):
            sys.stderr.write("{} does not exists. Please check input --counts\n".format(args.counts))
            exit(-1)
    elif args.barcodes is not None:
        if not os.path.isfile(args.barcodes):
            sys.stderr.write("File {} does not exist (--barcodes)\n".format(args.barcodes))
            exit(-1)
    else:
        sys.stderr.write("Please provide input --counts (mRNA data table) "
                         "or --barcodes (list of cell labels, no header).\n")
        exit(-1)

    if os.path.exists(args.output):
        sys.stderr.write("{} exists\n".format(args.output))
        exit(-1)
    os.mkdir(args.output)

    # count the tags, or use the counts provided
    if args.bam is not None:
        logging.info("Counting tags")
        cnt_df = count_tags_in_bam(args.bam, args.min_cell_sample_read, args.tag_gtf, cell_tag='CB')
    elif args.tag_counts is not None:
        cnt_df = bdgutils.mtx_split.GeneCellMatrix.factory(args.tag_counts).df.transpose()
        cnt_df.index.rename('cell', inplace=True)
    elif args.tag_ids is not None:
        pass;  # handle below after reading counts file
    else:
        sys.stderr.write("bam or tag_counts must be provided\n")
        exit(-1)

    # Get Results to be split and use them to filter the cells before
    # determining the sample assignments
    if args.counts is not None:
        gcm = bdgutils.mtx_split.GeneCellMatrix.factory(args.counts)
        filtered_barcodes = pd.DataFrame({'cell': gcm.cells})
    else:
        filtered_barcodes = pd.read_csv(args.barcodes, header = None, names = ['cell'])

    if args.tag_ids is not None:
        tag_ids = args.tag_ids.split(',')
        logging.info("tag_ids = %s", tag_ids)
        if args.counts is not None:
            cnt_df = gcm.gene_subset(tag_ids)
            cnt_df.index.rename('cell', inplace=True)
            logging.info("cnt_df = %s", cnt_df)

    logging.info("cells before filtering %d", len(cnt_df))
    logging.debug("cells : %s", cnt_df.index)
    logging.info("num filtered cells: %d", len(filtered_barcodes))
    logging.debug("filtered_cells: %s", filtered_barcodes)

    # error handling for number of cells found in bam file
    if len(cnt_df)==0:
        message = 'Could not find any reads mapped to sample tags. Please check that the correct gtf info is supplied.'
        raise Exception(message)

    if args.save_raw_counts:
        cnt_df.to_csv(os.path.join(args.output, "tagReadCounts.csv"))

    filtered_cntdf = pd.merge(filtered_barcodes, cnt_df, how='left',
                              left_on='cell', right_index=True).set_index('cell')

    assigner = dict()
    original_columns = filtered_cntdf.columns
    column_sums = dict()
    for col in original_columns:
        column_sums[col] = filtered_cntdf[col].sum()
    dropped_columns = [col for col in original_columns if column_sums[col] <= 1000]
    logging.info("Sum of sample tag counts: %s", column_sums)
    logging.debug("drop cols R = %s", dropped_columns)

    # check number of columns left
    if len(dropped_columns)==len(original_columns):
        message = 'Sample tag read counts are too low for analysis. Please refer to the log file {}/demux_counts.log for details.'.format(os.getcwd())
        raise Exception(message)

    # assign cells to a sample. Multiple methods can be used (specified by --method)
    logging.debug("original_columns = %s", original_columns)
    assigner[args.method] = SampleChooser.factory(args.method, filtered_cntdf[original_columns])
    filtered_cntdf[args.method] = assigner[args.method].assign_cells()

    filtered_cntdf.to_csv(os.path.join(args.output, "cell2Label.csv"))

    # get sample tag labels based on method
    sample_tag_calls = filtered_cntdf[[args.method]]
    sample_tag_calls.to_csv(os.path.join(args.output, "Sample_Tag_Calls.csv"))

    # generate and output the stats
    with open(os.path.join(args.output, 'stats.csv'), 'w') as outfile:
        write_summary_stats(outfile,
                            assigner,
                            *compute_stats(cnt_df, filtered_cntdf, assigner))


    # split the gene-cell matrix file by sample
    if args.counts is not None:
        logging.info("Splitting file to %s based on %s",
                     args.output, os.path.join(args.output, "cell2Label.csv"))
        gcm.set_sample_barcode_map(filtered_cntdf.reset_index().set_index(args.method))
        gcm.set_output_dir(args.output)
        gcm.write_per_sample_data()

    # output message
    print('Run completed. Output files are stored in {}'.format(args.output))


if __name__ == "__main__":
    sys.exit(main())
