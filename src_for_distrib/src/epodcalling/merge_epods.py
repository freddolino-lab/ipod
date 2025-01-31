#!/usr/bin/python

# calculate peak calls for all ipod samples of interest
# we read all of the instances to look at from a table of name/config file pairs

import argparse
import os
import sys
import numpy as np

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../utils")
sys.path.insert(0, utils_path)

import hdf_utils
import anno_tools as anno
import peak_utils as pu

def load_for_epod_merge(fnames):
    '''Reads in a list of narrowpeak files and returns a dictionary,
    the keys of which are indices of the file names in fnames,
    the values to which are anno.NarrowPeakData objects.

    Args:
    -----
    fnames : list
        List of narrowpeak file names to read in
    
    Returns:
    --------
    np_dict : dict
    '''
    
    np_dict = {}
    contigs = False
    for i,fname in enumerate(fnames):
        np_data = anno.NarrowPeakData()
        np_data.parse_narrowpeak_file(fname)
        if contigs:
            contigs.union(np_data.ctg_names())
        else:
            contigs = np_data.ctg_names()
        np_dict[i] = np_data

    return np_dict,contigs

def add_contig_records_to_np(contig_np, np_dict, ctg_id):
    for idx,np_data in np_dict.items():
        for record in np_data:
            if record.chrom_name == ctg_id:
                # set the record's name to idx so we can track its source later
                record.name = str(idx)
                contig_np.add_entry(record)
    contig_np.sort()
 
def merge_epods(fnames):
    '''Merges the regions found in the narrowpeak files listed in fnames.

    Args:
    -----
    fnames : list
        List of narrowpeak file names containing epods to merge.

    Returns:
    --------
    merged_epods : anno_tools.NarrowPeakData
        A NarrowPeakData object containing merged epods.
    '''

    np_dict,contigs = load_for_epod_merge(fnames)

    merged_epods = anno.NarrowPeakData()
    for ctg_id in contigs:

        ctg_np = anno.NarrowPeakData()
        # modify ctg_np in place here, adding records from the NarrowPeakData
        #  objects in the values of np_dict belonging to this ctg_id to ctg_np
        # ctg_np will then contain ALL records from ALL replicates from THIS
        #  contig. Replicate indices are tracked as the "name" attribute for
        #  each record.
        add_contig_records_to_np(ctg_np, np_dict, ctg_id)
 
        merged_contig_epods = merge_epods_in_contig(ctg_np)
        for epod in merged_contig_epods:
            merged_epods.add_entry(epod)

    return merged_epods


def get_grouped_intervals(narrowpeak):

    grpd_intervals = [[]]
    curr_start = narrowpeak[0].start
    curr_stop = narrowpeak[0].end
    samp_names = []

    for record in narrowpeak:
        if not record.name in samp_names:
            samp_names.append(record.name)
        if record.start < curr_stop:
            curr_stop = max(record.end, curr_stop)
            grpd_intervals[-1].append((record, record.name))
        else:
            curr_start = record.start
            curr_stop = record.end
            grpd_intervals.append([(record, record.name)])

    n_samples = len(samp_names)
    return grpd_intervals,n_samples

def merge_grouped_intervals(intervals, n_samples):

    # instantiate NarrowPeakEntry to store epod info
    epod = anno.NarrowPeakEntry()
    # start with absurd start/end values
    epod.start,epod.end = 1e12, -1
    # We'll keep track of the total length of all replicates' records in
    #  this merged epod
    num_reps_present = 0
    rep_inds = []
    rec_lengths = []
    rec_scores = []
    for record,rep_idx in intervals:

        if not rep_idx in rep_inds:
            rep_inds.append(rep_idx)
            num_reps_present += 1

        rec_lengths.append(record.end - record.start)
        rec_scores.append(record.score)
        epod.start = int(np.min([record.start, epod.start]))
        epod.end = int(np.max([record.end, epod.end]))

    epod.chrom_name = record.chrom_name
    epod.name = '.'

    # now we can calculate the average fraction of the epod represented
    #  by our replicates.
    epod_width = epod.end - epod.start
    record_cumulative_len = np.sum(rec_lengths)
    epod_replicate_representation = record_cumulative_len / epod_width
    epod_mean_rep_frac = epod_replicate_representation / n_samples
    # it is also usefult to know simply how many replicates had an
    #  epod in this merged epod
    epod_abs_rep_frac = num_reps_present / n_samples
    # Finally, we'll want to calculate the signal in the merged epod.
    #  We can do this by taking a weighted mean of each contributing
    #  epod's signal from each replicate. The weight is simply the
    #  fraction of the total merged epod length contained in this epod.
    # Get each record's fractional length
    rec_frac_len = np.asarray(rec_lengths) / epod_width

    # Now multiply each score by its fraction, and take the sum to get
    #  the weighted mean score for this merged epod.
    numer = np.sum(
        np.asarray(rec_scores)
        * rec_frac_len
    )
    denom = np.sum(rec_frac_len)
    mean_score = numer / denom
    
    epod.score = mean_score
    epod.qval = epod_mean_rep_frac
    epod.pval = epod_abs_rep_frac

    return epod
             
def build_unified_epods_narrowpeak(grouped_intervals, n_samples):

    merged_epods = []
    for intervals in grouped_intervals:
        merged_epod = merge_grouped_intervals(intervals, n_samples)
        merged_epods.append(merged_epod)

    return merged_epods

def merge_epods_in_contig(narrowpeak):
    '''Iterates over the narrow peak entires in narrowpeak, grouping
    overlapping intervals.
    '''

    # sort the entries first
    narrowpeak.sort()
    # grpd_intervals is a list of lists. The inner lists contain
    #  a tuple of (NarrowPeakEntry, replicate_idx) for each epod
    #  that overlaps in a single merged EPOD region.
    grpd_intervals,n_reps = get_grouped_intervals(narrowpeak)
    merged_epods = build_unified_epods_narrowpeak(grpd_intervals, n_reps)

    return merged_epods


if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--infiles',
        nargs = '+',
        help = "space-separated list of input narrowpeak files containing peaks or epods."
    )
    parser.add_argument(
        '--outfile',
        help = "Name of the narrowpeak file in which to place merged epods or peaks."
    )
    args = parser.parse_args()

    merged_epods = merge_epods(args.infiles)
    merged_epods.fname = args.outfile
    merged_epods.write_file()
