#!/usr/bin/python

import os
import sys
import numpy as np
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
    "--out_prefix",
    help = "Prefix, including path, to append final elements of output file names to.",
    required = True,
    type = str,
)
parser.add_argument(
    "--invert_scores",
    help="look for regions of very low occupancy instead of very high. Adds _LPOD to all filenames and acts on negative occupancy",
    action="store_true",
)

args = parser.parse_args()

OUTPREF = args.out_prefix

# read in bedgraph file
bg_info = anno.BEDGraphData()
bg_info.parse_bedgraph_file(args.in_file)
# sort by contig name and start position
bg_info.cleanup()
# get distinct contig names
ctgs = bg_info.ctg_names()

epod_results = anno.NarrowPeakData()

# loop over contigs to call peaks
for ctg_id in ctgs:

    ctg_info = anno.BEDGraphData()
    for record in bg_info:
        if record.filter('chrom_name', ctg_id):
            ctg_info.add_entry(record)
   
     
    
print("\n================================================")
print("Writing to {}".format(args.out_file))
results.write_file(args.out_file)
print("------------------------------------------------\n")

