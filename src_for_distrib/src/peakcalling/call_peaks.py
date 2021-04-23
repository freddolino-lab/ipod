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
import hdf_utils

parser = argparse.ArgumentParser()
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
#NOTE: change to --in_file
parser.add_argument(
    "--hdf_file",
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
###########################################################################
###########################################################################
###########################################################################
###########################################################################
#NOTE: deprecated
parser.add_argument(
    "--dataset_str",
    help = "Text appearing in the dataset name.",
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

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
#NOTE: get rid of all hdf5 file dependencies
# get contig lookup table
ctg_lut = hdf_utils.get_ctg_lut(args.hdf_file)
dset_fmt_str = args.dataset_str

# loop over contigs to call peaks
for ctg_id,ctg_info in ctg_lut.items():

    # load in data
    dset_name = dset_fmt_str.format(ctg_id)
    ctg_data = hdf_utils.load_dset(args.hdf_file, dset_name)

    kern = np.expand_dims(np.ones(args.window_size) / args.window_size, -1)

    rollmeans = scipy.signal.convolve(ctg_data, kern, mode="same")
    # label loci passing threshold as 1, others as 0
    goodflags = 1 * (rollmeans > args.threshold)
    num_good = np.sum(goodflags == 1)

    print("There are {} positions in contig {} greater than the current threshold of {}.".format(num_good, ctg_id, args.threshold))

    # write result to hdf5 and bedgraph outputs
    grp_name = '/'.join(dset_name.split('/')[:-1])
    out_dset_basename = "{}_cutoff_{}_peaks".format(
        dset_name.split('/')[-1],
        args.threshold,
    )

    hdf_utils.write_dset(
        args.hdf_file,
        out_dset_basename,
        goodflags,
        goodflags.dtype,
        grp_name,
    )
    
superctg_arr = hdf_utils.concatenate_contig_data(
    args.hdf_file,
    dset_basename = "{}/{}".format(args.sample_type, out_dset_basename),
)
print("================================================")
print("Writing to {}".format(args.out_file))
hdf_utils.write_bedgraph(superctg_arr, args.hdf_file, args.out_file)
print("------------------------------------------------\n")

