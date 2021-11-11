#!/usr/bin/python

import numpy as np
import scipy.integrate
import scipy.stats
import os
import sys
import toml

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, 'utils')
sys.path.insert(0, utils_path)

import hdf_utils

import matplotlib
# allow saving of fig without x forwarding
matplotlib.use("Agg")
from matplotlib import pyplot as plt

# calculate chip enrichment on each input file in my sequencing manifest
# we plot , and calculate the AUC for, the curve of normalized enrichment (against the maximum observed) vs cumulative distribution
# we assume that bootstrapping has already run

conf_file_global = sys.argv[1]
conf_dict_global = toml.load(conf_file_global)
HDF = sys.argv[2]
# still not used, but could be useful for compiling a bunch of qc into a single file.
all_qc_fname = sys.argv[3]

# useful constants
QCDIR = conf_dict_global["qc"]["qc_direc"]
# directory containing the bootstrap replicates and the actual coverage
BSDIR = conf_dict_global["bootstrap"]["bootstrap_direc"]
#q_suffix = conf_dict_global["bootstrap"]["orig_suffix"] + ".npy" # suffix for the file to be considered

def get_enr_stats(infile, outpng, outtxt):
    '''Reads an input numpy array containing coverage,
    calculates the coverage distribution over percentiles,
    writes an image and file with quantile-based
    scores and area under the curve.

    Args:
    -----
    infile : str
        Path to the numpy array of interest
    outpng : str
        Path to save output png.
    outtxt : str
        Path to save output auc information file.

    Returns:
    --------
    None
    '''

    # get dict of contig info from hdf file
    ctg_lut = hdf_utils.get_ctg_lut(infile)
    # grab first contig name from keys in ctg_lut
    this_ctg = [ctg_id for ctg_id in ctg_lut.keys()][0]
    group_name = "contigs/{}".format(this_ctg)

    # get values for just a single contig
    vals = hdf_utils.load_dset(infile, "orig", group_name)
    vals = vals/np.max(vals)

    quantiles = np.arange(100)

    scores = np.array( [
        scipy.stats.scoreatpercentile(vals, q)
        for q in quantiles
    ] )

    plt.plot(quantiles,scores)
    plt.savefig(outpng)

    # also calculate and write the auc
    auc = scipy.integrate.trapz(x=quantiles,y=scores)
    fout = open(outtxt,'w')
    fout.write("{}\n".format(auc))
    fout.close()


if __name__ == '__main__':

    if not os.path.isdir(QCDIR):
        os.mkdir(QCDIR)

    # NOTE: low priority, could come back to this to make the
    #   suffix a conf option
    this_prefix = os.path.basename(HDF).split('.')[0]
    out_img = os.path.join(QCDIR, this_prefix + '_chipqc.png')
    out_txt = os.path.join(QCDIR, this_prefix + '_chipqc.txt')
    get_enr_stats(HDF, out_img, out_txt)

