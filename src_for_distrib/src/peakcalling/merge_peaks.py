#!/usr/bin/python

'''Writes a file containing all regions of the given genome for which
at least the required fraction of input files have an overlapping region.
'''

import numpy as np
import argparse
import os
import sys

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../../src/utils")
sys.path.insert(0, utils_path)

import peak_utils as pu
import hdf_utils

if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--infiles',
        nargs = '+',
        help = "space-separated list of input narrowpeak files containing peaks.",
        required = True,
    )
    parser.add_argument(
        '--outfile',
        help = "Output narrowpeak file name.",
        required = True,
        type = str,
    )
    parser.add_argument(
        '--ref_db',
        help = "Absolute path to bowtie2 index of reference genome.",
        type = str,
        required = True,
    )
    parser.add_argument(
        '--threshold',
        help = "sets the threshold fraction of files in which a location must be represented to be included as a region in the output narrowpeak file.",
        required = True,
        type = float,
    )
    parser.add_argument(
        '--resolution',
        help = "The resolution at which analysis was performed.",
        required = True,
        type = int,
    )
    args = parser.parse_args()

    # get genome size from bowtie2 index
    ctg_lut = hdf_utils.make_ctg_lut_from_bowtie(args.ref_db)

    peaks_np = pu.merge_peaks_by_vote(
        args.infiles,
        ctg_lut,
        args.resolution,
        args.threshold,
    )
    peaks_np.fname = args.outfile
    peaks_np.write_file()

