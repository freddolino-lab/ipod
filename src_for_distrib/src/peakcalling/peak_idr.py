#!/usr/bin/python

import os
import sys

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import anno_tools as anno

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
    help = "Path to bedgraph file to write data to.",
    required = True,
    type = str,
)
parser.add_argument(
    "--window_size",
    help = "Size of the window to convolved rolling mean over",
    required = True,
    type = int,
)
parser.add_argument(
    "--threshold",
    help = "Threshold value at which to call a locus a peak\
            after taking rolling mean.",
    required = True,
    type = float,
)

args = parser.parse_args()

# read in bedgraph file
bg_info = anno.BEDGraphData()
bg_info.parse_bedgraph_file(args.in_file)
# sort by contig name and start position
bg_info.cleanup()
# get distinct contig names
ctgs = bg_info.ctg_names()

results = anno.NarrowPeakData()

# loop over contigs to call peaks
for ctg_id in ctgs:

    ctg_info = anno.BEDGraphData()
    for record in bg_info:
        if record.filter('chrom_name', ctg_id):
            ctg_info.add_entry(record)
    
    # get the score for each position
    scores = np.expand_dims(ctg_info.fetch_array(attr='score'), -1)
    starts = ctg_info.fetch_array(attr='start')
    ends = ctg_info.fetch_array(attr='end')

    kern = np.expand_dims(np.ones(args.window_size) / args.window_size, -1)

    rollmeans = scipy.signal.convolve(scores, kern, mode="same")
    # label loci passing threshold as 1, others as 0
    goodflags = 1 * (rollmeans > args.threshold)

    in_peak = False
    state_change = False
    peak_scores = []

    for i,flag in enumerate(goodflags):

        # decide whether we switched state
        if not in_peak:
            if flag == 1:
                state_change = True
        else:
            if flag == 0:
                state_change = True

        # determine if we're in a peak
        if flag == 1:
            in_peak = True
        else:
            in_peak = False

        # if we moved out of a peak, record prior end position and calculate peak
        #   score and point-source
        if (not in_peak) and state_change:
            # index is for prior site, since this one is just outside the peak.
            #   Subtract 1 since narrowpeak is zero-indexed
            peak_end = ends[i-1]
            peak_score = np.mean(peak_scores)
            ####################################################
            ####################################################
            ######### calculate peak point-source too
            results.addline(
                chrom_name = ctg_id,
                start = peak_start,
                end = peak_end,
                score = peak_score,
            )
            
            peak_scores = []
            state_change = False

        # if we are in a peak, record score
        if in_peak:
            peak_scores.append(rollmeans[i])
            
            # if we moved into a peak, record start and scores within peak
            if state_change:
                peak_start = starts[i]
                state_change = False

    num_peaks = len(results)

    print("There are {} peaks in contig {} passing the current threshold of {}.".format(num_peaks, ctg_id, args.threshold))
    
print("\n================================================")
print("Writing to {}".format(args.out_file))
results.write_file(args.out_file)
print("------------------------------------------------\n")

