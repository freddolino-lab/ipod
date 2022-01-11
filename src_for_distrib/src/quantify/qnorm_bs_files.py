#!/usr/bin/python

import h5py
import numpy as np
import scipy.stats
import toml
import os
import sys
import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import quant_utils as qutils
import hdf_utils

# Run quantile normalization on all bootstrap output files.
# We do the quantile normalization separately for the each signal type
#   defined in our config file.
# note that all distributions - both the original and bootstrap - are
#   normalized based on the quantiles from the median of the original
#   distributions for that sample type.

def calc_qnorm_base( input_arrays ):
    '''Calculate the averaged quantile-wise distribution 
    over a bunch of targets. We take the median of the input 
    values at each quantile to define that quantile.
    Note that we assume that there are no missing values, 
    and that all inputs have the same length.

    Args:
    -----
    input_arrays : list
        A list containing each array of coverages

    Returns:
    --------
    median_qs : numpy array
        A 1d numpy array containing the median coverage
        for each position among all samples in the input_arrays
    '''


    # first we sort the inputs column-wise and place sorted vals into
    # an array
    sorted_arrs = np.zeros((input_arrays[0].shape[0], len(input_arrays)))
    for i in range(len(input_arrays)):
        this_arr = input_arrays[i]
        sort_inds = np.argsort(this_arr[:,0], kind='heapsort')
        sorted_arrs[:,i] = this_arr[sort_inds,0]

    # now get the median across samples at each position
    median_qs = np.median(sorted_arrs, axis=1)

    return median_qs

def q_norm_vec( input_vals, target_vals ):
    '''Quantile normalize input_vals to match the 
    distribution in target_vals.

    Args:
    -----
    input_vals : numpy array
        1d numpy array containing coverage at each position
    target_vals : numpy array
        
    Returns:
    --------
    rank-ordered target_vals : numpy array
        this array contains the original target vals, ranked according to
        the sorted order of the original data in the input_vals array.
        This means that the rank-ordered target vals are the quantile
        normalized data for this sample.
    '''

    rank_vec = scipy.stats.rankdata(input_vals, method='ordinal')
    qnorm_vals = target_vals[rank_vec - 1]
    # divide lowest non-zero value by div, and add the result
    #  to each zero value.
    qnorm_vals = qutils.add_pseudocount_vec(
        vec = qnorm_vals,
        div = 2,
    )

    return qnorm_vals

def qnorm_bootstrap_mat(full_mat, target_dist):
    '''Apply bootstrap normalization separately to each column in
    a matrix, using target_dist as the target set of values.
    If outfile is None, overwrite the input.

    Args:
    -----
    full_mat : 2d np.array
        Array of shape (genome_length,nboot).
    target_dist : 1d np array
        Array containing the target values for quantile normalization.

    Modifies:
    ---------
    full_mat
        Modified in place.
    '''

    nvals,nboot = full_mat.shape

    for b in range(nboot):
        this_vec = full_mat[:,b]
        full_mat[:,b] = q_norm_vec(this_vec, target_dist)

def q_norm_files(hdf_names, ctg_lut, out_dset_name, bs_num):
    '''Read and quantile normalize a set of files. We read all of the 
    members of orig_files, calculate the target distribution based
    on them, and then write the normalized version of each. Then, 
    every bootstrap replicate in each of the files in bs_files is normalized 
    TO THE SAME TARGET DISTRIBUTION.

    Args:
    -----
    hdf_names : list
        Names of files containing actual coverage.
    ctg_lut : dict
        Dictionary containing information on which contig in the reference
        sequence corresponds to information in the hdf file.
    out_dset_name : str
        Name of new dataset to create in each hdf file.
    bs_num : int
        Number of bootstrap samples taken at bootstrap time.

    Returns:
    --------
    None
        Writes output as new dataset in hdf5 files.
    ''' 

    # Calculate the target distribution and rewrite the original files.
    orig_vecs = []
    # loop over each sample's data, appending each contig's data into
    #   one long supercontig, and append coverage to list.
    for fname in hdf_names:
        print(fname)
        concat_arr = hdf_utils.concatenate_contig_data(
            fname,
        )
        orig_vecs.append(concat_arr)

    print('calculating target distribution.......')
    # target_distr has the medians
    target_distr = calc_qnorm_base(orig_vecs)

    print('quantile normalizing original data......')
    qnorm_vecs = [
        np.expand_dims(q_norm_vec( v, target_distr), -1)
        for v in orig_vecs
    ]

    for i,fname in enumerate(hdf_names):

        hdf_utils.decatenate_and_write_supercontig_data(
            fname,
            qnorm_vecs[i],
            dset_name = "orig_{}".format(out_dset_name),
            attrs = {"normalization_method": "quantile"},
        )

    # now do the same for each of the bootstrap files
    print('quantile normalizing bootstrap data.....')
    for fname in hdf_names:

        these_vals = hdf_utils.concatenate_contig_data(
            fname,
            dset_basename = "bs",
            sample_num = bs_num,
        )

        # these_vals modified in place here.
        qnorm_bootstrap_mat(these_vals, target_distr)

        hdf_utils.decatenate_and_write_supercontig_data(
            fname,
            these_vals,
            dset_name = "bs_{}".format(out_dset_name),
            dtype = np.float64,
            attrs = {"normalization_method": "quantile"},
        )
        
# here is where the main program starts

if __name__ == '__main__':

    conf_file = sys.argv[1]
    conf_dict = toml.load(conf_file)
    conf_file_global = sys.argv[2]
    conf_dict_global = toml.load(conf_file_global)

    # figure out some global parameters
    BS_DIR = conf_dict_global['bootstrap']['bootstrap_direc']
    BS_NUM = conf_dict_global['bootstrap']['bootstrap_samples']
    OUT_DSET = conf_dict_global['norm']['qnorm_dset']
    out_prefix = os.path.join(
        conf_dict_global['bootstrap']['output_path'],
        conf_dict['general']['out_prefix']
    )

    # run quantile normalization on each set of samples defined in
    #  the config file
    sample_types = conf_dict["general"]["sample_types"]
    qnorm_samples = conf_dict["quant"]["qnorm_samples"]

    for samptype in sample_types:
        # skip if this sample type is not in the qnorm samples list
        if not samptype in qnorm_samples:
            continue
        data_dir = conf_dict[samptype]['directory']
        sample_prefixes = conf_dict[samptype]['sample_names']
        #print(sample_prefixes)
        sample_hdf_fnames = [
            os.path.join( data_dir, BS_DIR, pref + '.hdf5' )
            for pref in sample_prefixes
        ]

        # Set up ctg_lut to organize ctg ids and indices in arrays
        ctg_lut = hdf_utils.get_ctg_lut(sample_hdf_fnames[0])
        q_norm_files(
            sample_hdf_fnames,
            ctg_lut,
            OUT_DSET,
            BS_NUM,
        )

