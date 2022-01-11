#!/usr/bin/env python3

"""
Given a set of narrowpeak files containing peaks identified at
the given thresholds, choose the best threshold for
discriminating between peak and non-peak regions.
"""

import argparse
import shutil
import os
import sys
import toml
import subprocess
import glob
import re
import numpy as np
import tempfile
import multiprocessing
from recombinator import iid_bootstrap
import scipy.stats
from matplotlib import pyplot as plt

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../utils")
sys.path.insert(0, utils_path)

import hdf_utils
import anno_tools as anno
import peak_utils as pu

def get_kl_divergences(peak_scores, nonpeak_scores,
                    lowest_bin=-10, highest_bin=11, alpha=0.05, n_b=1000):
    """Calculate KL divergence between peaks and non-peaks. Also
    performs bootstrapping to get the upper and lower confidence limits
    for the KL divergence.
    """

    if peak_scores.size == 0:
        return (0,0,0)
    elif peak_scores.size == 1:
        peak_scores = np.expand_dims(peak_scores, -1)

    bins = np.arange(lowest_bin, highest_bin)
    counts_1,bins_1 = np.histogram(peak_scores, bins=bins)
    counts_2,bins_2 = np.histogram(nonpeak_scores, bins=bins)

    entropy=scipy.stats.entropy( counts_1+1, counts_2+1 )

    value_distr_1 = iid_bootstrap(peak_scores, replications=n_b, replace=True)
    value_distr_2 = iid_bootstrap(nonpeak_scores, replications=n_b, replace=True)

    #print(len(peak_scores))
    #print(n_b)
    #value_distr_1 = np.zeros((n_b, len(peak_scores)))
    #value_distr_2 = np.zeros((n_b, len(nonpeak_scores)))

    ## generate bootstraps
    #for r in range(value_distr_1.shape[0]):
    #    value_distr_1[r,:] = random.choices(peak_scores, k=n_b)
    #    value_distr_2[r,:] = random.choices(nonpeak_scores, k=n_b)

    entropy_from_bootstrap=[]

    for i in range(n_b):
        counts_1, bins_1 = np.histogram(
            value_distr_1[i,:],
            bins=bins,
        )
        counts_2, bins_2 = np.histogram(
            value_distr_2[i,:],
            bins=bins,
        )
        entropy_boot = scipy.stats.entropy( counts_1+1, counts_2+1 )
        entropy_from_bootstrap.append( entropy_boot )

    # now get the confidence interval based on percentiles
    cl_lo = scipy.stats.scoreatpercentile( entropy_from_bootstrap, 100*(alpha/2) )
    cl_hi = scipy.stats.scoreatpercentile( entropy_from_bootstrap, 100*(1- (alpha/2)) )

    return (cl_lo, entropy, cl_hi)

           
def choose_final_threshold(idr_files, ctg_lut, spike_name, mean_fname,
                        lowest_bin=-10, highest_bin=11, alpha=0.05, n_b=1000,
                        debug=False):
    """Get the cutoff index that we consider to be the best cutoff for
    robustly distinguishing between peak and non-peak regions of the genome.
    """

    divergences = []
    #print(idr_files)
    for i,idr_fname in enumerate(idr_files):
        print("=============================")
        print("Calculating KL divergence between peaks in {} and non-peaks.".format(idr_fname))
        print("=============================")

        final_peaks = anno.NarrowPeakData()
        final_peaks.parse_narrowpeak_file(idr_fname)
        complement_bed_data = final_peaks.fetch_complement_bed_data(
            contig_lut = ctg_lut,
            filter_chrs = [spike_name],
        )

        # write the bed file containing peak-less regions to bed file
        if debug:
            nonpeak_bed_fname = '/home/schroedj/nopeak_data.bed'
            complement_bed_data.fname = nonpeak_bed_fname
            complement_bed_data.write_file()
            peak_score_fname = '/home/schroedj/peakscore.bed'
            nonpeak_score_fname = '/home/schroedj/nonpeakscore.bed'

            print("Mean file name: {}".format(mean_fname))
            print("Peak score file name: {}".format(peak_score_fname))
            print("Non-peak score file name: {}".format(nonpeak_score_fname))
            print("Non-peak regions file name: {}".format(complement_bed_data.fname))
        else:
            tmp_dir = tempfile.TemporaryDirectory()
            nonpeak_bed_fname = os.path.join(tmp_dir.name, 'nopeak_data.bed')
            complement_bed_data.fname = nonpeak_bed_fname
            complement_bed_data.write_file()

            # write mean peak and non-peak scores to temporary files
            peak_score_fname = os.path.join(tmp_dir.name, 'peakscore.bed')
            nonpeak_score_fname = os.path.join(tmp_dir.name, 'nonpeakscore.bed')

        bed_map_cmd = "bedtools map \
                -a {} \
                -b {} \
                -o mean -c 4 \
                | cut -f 7 \
                > {}".format(
                final_peaks.fname,
                mean_fname,
                peak_score_fname,
        )
        subprocess.call(
            bed_map_cmd,
            shell = True,
        )

        peak_scores = np.loadtxt(peak_score_fname)
        
        bed_map_cmd = "bedtools map \
                -a {} \
                -b {} \
                -o mean -c 4 \
                | cut -f 7 \
                > {}".format(
                nonpeak_bed_fname,
                mean_fname,
                nonpeak_score_fname,
        )
        subprocess.call(
            bed_map_cmd,
            shell = True,
        )

        nonpeak_scores = np.loadtxt(nonpeak_score_fname)

        if not debug:
            tmp_dir.cleanup()

        # give a tuple of (lower_cl, observed, upper_cl) for the KL divergence
        this_div_info = get_kl_divergences(
            peak_scores,
            nonpeak_scores,
            lowest_bin,
            highest_bin,
            alpha,
            n_b,
        )
        divergences.append(this_div_info)
    # make divergences an array, and transpose it and slice rows in reverse
    #  so final array's rows are upper, observed, lower at indices 0,1,2, repectively,
    #  and columns are cutoffs
    try:
        div_arr = np.asarray(divergences).T[::-1,:]
    except:
        print(divergences)
        sys.exit()
    max_observed = div_arr[1,:].max()
    # get max index at which upper cl is greater than the max observed KL divergence
    cutoff_idx = np.where(div_arr[0,:] > max_observed)[0].max()
    return (cutoff_idx, div_arr, max_observed)
 
def plot_results(best_thresh_path, basename, kl_divs, cutoffs, size=(9,5)):

    plt_name = os.path.join(
        best_thresh_path,
        "kl_div_cutoff_{}.png".format(basename),
    )

    # get error bar info as ax.errorbar requires (distance from observed)
    err_bar_arr = kl_divs[[2,0],:] - kl_divs[1,:][None,:]
    err_bar_arr[0,:] *= -1

    fig, ax = plt.subplots(figsize=size)
    # plot horizontal dashed red line at max observed kl_div
    ax.axhline(
        y=kl_divs[1,:].max(),
        color='r',
        linestyle='--',
    )
    # draw line with conf limits as bars
    ax.errorbar(
        x = cutoffs,
        y = kl_divs[1,:],
        yerr = err_bar_arr,
        color='tab:blue',
        ecolor='tab:blue',
    )
    ax.set_xlabel('threshold')
    ax.set_ylabel('KL divergence')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(plt_name)
    plt.close()


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--infiles',
        nargs = '+',
        help = "space-separated list of input narrowpeak files containing peaks.",
        required = True,
    )
    parser.add_argument(
        '--thresholds',
        nargs = '+',
        help = "space-separated list of thresholds used to call peaks in each input narrowpeak file.",
        required = True,
    )
    parser.add_argument(
        '--out_dir',
        help = "Directory into which to write output files.",
        type = str,
        required = True,
    )
    parser.add_argument(
        '--lowest',
        help = "Used for binning of scores in peak and non-peak regions. Sets the lowest bin value.",
        type = float,
        required = True,
    )
    parser.add_argument(
        '--highest',
        help = "Used for binning of scores in peak and non-peak regions. Sets the highest bin value.",
        type = float,
        required = True,
    )
    parser.add_argument(
        '--ref_db',
        help = "Absolute path to bowtie2 index of reference genome.",
        type = str,
        required = True,
    )
    parser.add_argument(
        '--mean_score_fname',
        help = "Absolute path to bedgraph file containing the score to map to peaks and non-peaks.",
        type = str,
        required = True,
    )
    parser.add_argument(
        '--spikein_name',
        help = "Name of contig in reference used for spike-in normalization.",
        type = str,
        default = "None",
    )
    parser.add_argument(
        '--debug',
        help = "set this argument to use debugging mode.",
        action = "store_true",
    )
    parser.add_argument(
        '--sample_num',
        help = "sets the number of bootstrap samples. Default is 1000.",
        default = 1000,
        type = int,
    )
    parser.add_argument(
        '--alpha',
        help = "sets the confidence limit to campare, for each threshold, to the max observed KL divergence. Default is 0.95.",
        default = 0.95,
        type = float,
    )
    args = parser.parse_args()

    ctg_lut = hdf_utils.make_ctg_lut_from_bowtie(args.ref_db)

    best_cutoff_idx,kl_divs,max_obs_kl_div = choose_final_threshold(
        args.infiles,
        ctg_lut,
        args.spikein_name,
        args.mean_score_fname,
        args.lowest,
        args.highest,
        args.alpha,
        args.sample_num,
        args.debug,
    )

    best_cutoff = args.thresholds[best_cutoff_idx]
    best_result = args.infiles[best_cutoff_idx]
    basename = os.path.basename(best_result)

    plot_results(args.out_dir, basename, kl_divs, args.thresholds)

    shutil.copy(best_result, os.path.join(args.out_dir, basename))


if __name__ == '__main__':
    main()


