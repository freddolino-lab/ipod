#!/usr/bin/python

import os
import sys
import numpy as np
import operator
import scipy.signal
import argparse

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import anno_tools as anno
import peak_utils as pu

parser = argparse.ArgumentParser()

parser.add_argument(
    "--in_file",
    help = "Name of hdf5 file to read and write data from/to.",
    required = True,
    type = str,
)
parser.add_argument(
    "--sample_type",
    help = "Name of the sample (ipod, nc, etc.).",
    required = True,
    type = str,
)
parser.add_argument(
    "--out_file",
    help = "Path to narrowpeak file to write data to.",
    required = True,
    type = str,
)
parser.add_argument(
    "--window_size",
    help = "Size of the window (in base-pairs) to convolve rolling median over",
    required = True,
    type = int,
)
parser.add_argument(
    "--threshold",
    help = "Threshold value at which to call a locus a peak\
            after taking rolling median.",
    required = True,
    type = float,
)

args = parser.parse_args()

# read in bedgraph file
bg_info = anno.BEDGraphData()
bg_info.parse_bedgraph_file(args.in_file)
# get distinct contig names
ctgs = bg_info.ctg_names()

results = anno.NarrowPeakData()

# loop over contigs to call peaks
for ctg_id in ctgs:

    # keep just this contig's records
    ctg_info = anno.BEDGraphData()
    for record in bg_info:
        if record.filter('chrom_name', ctg_id):
            ctg_info.add_entry(record)
    
    # get the score for each position
    scores = ctg_info.fetch_array(attr='score')
    #scores = np.expand_dims(ctg_info.fetch_array(attr='score'), -1)
    starts = ctg_info.fetch_array(attr='start')
    ends = ctg_info.fetch_array(attr='end')

    resolution = ends[0] - starts[0]

    if args.window_size == resolution:
        rollmedians = scores
    else:
        rollmedians = pu.calc_ctg_running_median(
            scores,
            args.window_size,
            int(ends[0]-starts[0]),
            units_bp = True,
        )
    # label loci passing threshold as 1, others as 0
    goodflags = 1 * (rollmedians > args.threshold)

    # results updated in-place here, so nothing returned
    pu.get_peaks_from_binary_array(
        ctg_id,
        starts,
        ends,
        goodflags,
        results,
        rollmedians,
    )

    num_peaks = len(results.filter("chrom_name", operator.eq, ctg_id))

    print(
        f"There are {num_peaks} peaks in contig {ctg_id} "\
        f"passing the current threshold of {args.threshold}."
    )
    
print("\n================================================")
num_peaks = len(results)
if num_peaks > 0:
    print("Writing to {}".format(args.out_file))
    results.fname = args.out_file
    results.write_file()
else:
    print("No peaks, not writing any output.")
print("------------------------------------------------\n")

