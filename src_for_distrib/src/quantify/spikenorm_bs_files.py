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

import hdf_utils
import quant_utils as qutils

# Run spike-in normalization on all bootstrap output files.
def trimmed_mean(in_arr, upper_quant=0.8, lower_quant=0.2):
    '''Convenience function for getting trimmed mean of each
    column of in_arr.
    '''
    # Get lower and upper quintile thresholds for doing trimmed mean
    quants = np.quantile(in_arr, [lower_quant, upper_quant], axis=0)
    # which indices are between the chosen quantiles?
    keepers = np.logical_and(
        in_arr > quants[0,:],
        in_arr < quants[1,:],
    )
    # keep those that were between the chosen quantiles
    kept = []
    for i in range(in_arr.shape[1]):
        these_keepers = keepers[:,i]
        these_kept = in_arr[these_keepers, i]
        kept.append(these_kept)
    trimmed_arr = np.stack(kept, axis=1)
    # this number represents the mean number of reads
    #   allocated to spike-in
    trimmed_means = trimmed_arr.mean(axis=0)
    return trimmed_means


def clipped_mean(in_arr, res, clip_len=50):
    '''Convenience function for clipping ends of in_arr
    and calculating mean for each column of in_arr.
    '''

    # identify how many positions will be clipped
    clip_pos = clip_len // res
    if clip_pos == 0:
        arr = in_arr
    else:
        # slice everything between the clipped positions
        arr = in_arr[clip_pos:-clip_pos,:]
    # this number represents the mean number of reads
    #   allocated to spike-in for each replicate/strand
    # shape (R,S), where R is replicate, strand is S
    means = arr.mean(axis=0)
    return means


def spike_normalize(genome_counts_arr, spike_counts_arr,
                    mean_spike_arr, spike_amount, sample_cfu):
    '''Determine the amount of material per cfu represented by the coverage at
    each position (row) in genome_counts_arr.

    Args:
    -----
    genome_counts_arr : np.array
        2d numpy array of shape (P,T,S), where P is the number of genome positions
        considered in this analysis at the chosen resolution and S is the
        number of bootstrap samples (will be 1 if it's just raw coverage).
        Values in this array are simply read counts piling up at each position.
    spike_counts_arr : np.array
        2d numpy array of shape (K,B), where K is the number of spike-in positions
        and B is the number of bootstrap samples (will be 1 if it's just observed
        coverage).
    mean_spike_arr : np.array
        Numpy array of shape (T,S) containing the trimmed mean number of reads
        aligning to the positions of the spike-in "chromosome" for each
        sample type, T and strand, S.
    spike_amount : float
        How much spike-in was provided to each
        sample. Bear in mind whether you provided a biological or an
        in vitro spike-in and track your units accordingly.
    sample_cfu : list
        List of length T identifying the total colony forming units that
        went into preparing each sample.

    Returns:
    --------
    amount_per_cfu : np.array
        Numpy array of shape (P,T,S), the values of which are the amount
        of material per cfu represented by the sequencing coverage
        at each position, P, in each sample type, T, on each strand, S.
    '''

    spike_amount_per_read = (
        np.expand_dims(np.array(spike_amount), axis=-1)
        / mean_spike_arr
    )

    amount_per_cfu = (
        genome_counts_arr
        * np.expand_dims(spike_amount_per_read, axis=0)
        / np.expand_dims(sample_cfu, axis=-1)
    )

    return amount_per_cfu


def spike_norm_files(hdf_names, ctg_lut, out_dset_name, bs_num,
                     spike_chr, sample_cfu, sample_spikein_amount,
                     resolution, clip_length_bp, diagnostic_file_names,
                     orig_dset="orig"):
    '''Read and spike-in normalize a set of files. We read all of the 
    members of orig_files, calculate the target distribution based
    on them, and then write the normalized version of each. Then, 
    every bootstrap replicate in each of the files in bs_files is normalized 
    to its own bootsptrap-sampled spike in.

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
    spike_chr : str
        Name of the "chromosome" to which reads from spike-in
        should align. Default value is None.
    sample_cfu : list
        Total colony forming units for each sample.
    sample_spikein_amount : list
        List containing the nanograms of spikein added to each sample.
    resolution : int
        Resolution for this experiment
    clip_length_bp : int
        Number of base pairs to clip from the ends of the spike-in
        to calculate mean coverage over the spike-in.
    diagnostic_file_names : str
        Names of files to write fraction of alignments aligning to the genome
        and fraction of alignments aligning to spike-in.
    orig_dset : str
        Sets base dataset name to read from hdf5 file. Default is "orig".

    Returns:
    --------
    None
        Writes output as new dataset in hdf5 files.
    ''' 

    # Calculate the target distribution and rewrite the original files.
    orig_vecs = []
    spike_vecs = []

    # loop over each sample's data, appending each contig's data into
    #   one long supercontig, and append coverage to list.
    print("spike-in normalizing empirical coverage data.....")
    for i,fname in enumerate(hdf_names):

        all_ctg_arr = hdf_utils.concatenate_contig_data(
            fname,
            dset_basename = orig_dset,
        )
        orig_vecs.append(all_ctg_arr)

        if spike_chr is not None:
            spike_arr = hdf_utils.get_contig_data(
                fname,
                contig_name = spike_chr,
                dset_basename = orig_dset,
            )
            spike_vecs.append(spike_arr)

        total = all_ctg_arr.sum()

        ctg_lut = hdf_utils.get_ctg_lut(fname)
        diagnostic_file_str = "Fraction aligning to {}"
        diagnostic_file_fields = []
        diagnostic_file_data = []
        # iterate over all contigs
        for ctg_id,ctg_info in ctg_lut.items():
            # get summed coverage for this contig
            this_ctg_cov = np.sum(
                all_ctg_arr[ctg_info['start_idx']:ctg_info['end_idx']]
            )
            # append this contig to header fields
            diagnostic_file_fields.append(diagnostic_file_str.format(ctg_id))
            # append this contig's fractional allocation of coverage to data fields
            diagnostic_file_data.append( str(this_ctg_cov / total) )

        #spikein_sum = spike_arr.sum()

        #frac_spikein = spikein_sum / total
        #frac_genome = 1 - frac_spikein

        with open(diagnostic_file_names[i], 'w') as outf:
            outf.write(f"{','.join(diagnostic_file_fields)}\n")
            outf.write(f"{','.join(diagnostic_file_data)}\n")

    # stack then slice so they're now of shape (G,R,S),
    # where G is genome len, R is rep num, and S is strand num
    genome_count_arr = np.stack(orig_vecs, axis=1)[:,:,0,:]
    spikein_count_arr = np.stack(spike_vecs, axis=1)[:,:,0,:]

    clipped_means = clipped_mean(
        spikein_count_arr,
        resolution,
        clip_length_bp,
    )

    # divide the amount of spike-in by number of reads and cfus
    #  to get amount of material per read per cfu
    amount_per_cfu = spike_normalize(
        genome_count_arr,
        spikein_count_arr,
        clipped_means,
        sample_spikein_amount,
        sample_cfu,
    )

    for i,fname in enumerate(hdf_names):

        hdf_utils.decatenate_and_write_supercontig_data(
            fname,
            np.expand_dims(amount_per_cfu[:,i,:], axis=1),
            dset_name = f"orig_{out_dset_name}",
            attrs = {'normalization_method': "spike-in"},
        )

    # now do the same for each of the bootstrap files
    print('spike-in normalizing bootstrap data.....')
    for i,fname in enumerate(hdf_names):

        these_genome = hdf_utils.concatenate_contig_data(
            fname,
            dset_basename = "bs",
            sample_num = bs_num,
        )

        these_spikes = hdf_utils.get_spikein_data(
            fname,
            spikein_name = spike_chr,
            dset_basename = "bs",
        )
        bs_clipped_means = clipped_mean(
            these_spikes,
            resolution,
            clip_length_bp,
        )

        # because here we're broadcasting our math on a
        #   single sample's bootstrap replicates, just use
        #   that sample's spike-in amount and cfu count
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
        bs_amount_per_cfu = spike_normalize(
            these_genome,
            these_spikes,
            bs_clipped_means,
            sample_spikein_amount,
            sample_cfu[i],
        )

        bs_mean_per_cfu = np.mean(bs_amount_per_cfu, axis=1)
        # divide lowest non-zero value by div, and add the result
        #  to each zero value.
        bs_mean_per_cfu = qutils.add_pseudocount_vec(
            vec = bs_mean_per_cfu,
            div = 2,
        )

        # write data to hdf5 file
        hdf_utils.decatenate_and_write_supercontig_data(
            fname,
            bs_amount_per_cfu,
            dset_name = f"bs_{out_dset_name}",
            dtype = np.float64,
            attrs = {'normalization_method': "spike-in"},
        )
        # write data to hdf5 file
        hdf_utils.decatenate_and_write_supercontig_data(
            fname,
            bs_mean_per_cfu,
            dset_name = f"bs_mean_{out_dset_name}",
            dtype = np.float64,
            attrs = {'normalization_method': "spike-in"},
        )
        f_path,base_name = os.path.split(fname)
        base_name = base_name.split('.')[0] + "_spikenorm_mean.bedgraph"
        bg_fname = os.path.join(f_path, base_name)

        hdf_utils.write_bedgraph(
            bs_mean_per_cfu,
            fname,
            bg_fname,
        )

        total = these_genome.sum(axis=0)
        spikein_sum = these_spikes.sum(axis=0)

        frac_spikein = spikein_sum / total
        frac_genome = 1 - frac_spikein

        bs_diagnostic_fnames = []
        for diag_fname in diagnostic_file_names:
            direc,basename = os.path.split(diag_fname)
            basename = "bs_" + basename
            bs_diagnostic_fnames.append(os.path.join(direc,basename))
            
        print(f"Spike-in normalizing boostrap samples for {fname}")
        ctg_lut = hdf_utils.get_ctg_lut(fname)
        bs_diagnostic_file_fields = []
        bs_diagnostic_file_data = []
        # iterate over all contigs
        for ctg_id,ctg_info in ctg_lut.items():
            # get summed coverage for this contig, for each bootstrap
            this_ctg_cov = np.sum(
                these_genome[ctg_info['start_idx']:ctg_info['end_idx'],:],
                axis = 0,
            )
            # append this contig to header fields
            bs_diagnostic_file_fields.append( diagnostic_file_str.format(ctg_id) )
            # append this contig's fractional allocation of coverage to data fields
            bs_diagnostic_file_data.append( this_ctg_cov / total )

        with open(bs_diagnostic_fnames[i], 'w') as outf:
            outf.write(f"Bootstrap replicate,{','.join(bs_diagnostic_file_fields)}\n")
            for j in range(these_genome.shape[1]):
                this_bs_data = [str(ctg_vals[j]) for ctg_vals in bs_diagnostic_file_data]
                outf.write(f"{j+1},{','.join(this_bs_data)}\n")
        

# here is where the main program starts
if __name__ == '__main__':

    conf_file = sys.argv[1]
    conf_dict = toml.load(conf_file)
    conf_file_global = sys.argv[2]
    conf_dict_global = toml.load(conf_file_global)

    # figure out some global parameters
    STRANDED = conf_dict_global['general']['stranded']
    BS_DIR = conf_dict_global['bootstrap']['bootstrap_direc']
    BS_NUM = conf_dict_global['bootstrap']['bootstrap_samples']
    OUT_DSET = conf_dict_global['norm']['spikenorm_dset']
    RES = conf_dict_global['genome']['resolution']
    CLIP = conf_dict_global['norm']['clip_len_bp']
    SPIKE_CHR = conf_dict_global['genome']['spike_in_name']
    out_prefix = os.path.join(
        conf_dict_global['bootstrap']['output_path'],
        conf_dict['general']['out_prefix']
    )
    sample_spikein_amount = conf_dict['quant']['spikein_amount']
    cfu_scaling = conf_dict_global['general']['cfu_scale_fac']
    sample_cfu = np.asarray(conf_dict['quant']['cfu']) / cfu_scaling
    spikenorm_samples = conf_dict['quant']['spikenorm_samples']

    # run quantile normalization on each set of samples defined in
    #  the config file
    sample_types = conf_dict["general"]["sample_types"]

    for samptype in sample_types:
        if not samptype in spikenorm_samples:
            continue
        if SPIKE_CHR == "None":
            sys.exit("Error: you have not named your spike-in chromosome that is in your reference genome fasta file. You must provide its *exact* name as the ['genome']['spike_in_name'] option to your main config file located at {} to run spike-in normalization".format(conf_file_global))
 
        data_dir = conf_dict[samptype]['directory']
        sample_prefixes = conf_dict[samptype]['sample_names']

        sample_hdf_fnames = [
            os.path.join( data_dir, BS_DIR, pref + '.hdf5' )
            for pref in sample_prefixes
        ]
        diagnostic_fnames = [
            os.path.join( data_dir, BS_DIR, pref + "_read_allocation.csv" )
            for pref in sample_prefixes
        ]

        # Set up ctg_lut to organize ctg ids and indices in arrays
        ctg_lut = hdf_utils.get_ctg_lut(sample_hdf_fnames[0])

        spike_norm_files(
            sample_hdf_fnames,
            ctg_lut,
            OUT_DSET,
            BS_NUM,
            SPIKE_CHR,
            sample_cfu,
            sample_spikein_amount,
            RES,
            CLIP,
            diagnostic_fnames,
        )

