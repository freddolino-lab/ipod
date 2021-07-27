#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")

import numpy as np
from functools import partial
import toml
import os
import subprocess
import sys
import re
from pprint import pprint
import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import hdf_utils
import quant_utils as qutils

# Functions to apply ipod data processing and jackknife resampling.
# In the end this should allow us to collapse all of the replicates
#   from ipod, chip, and inp samples to a set of summary files:
#  -ipod_vs_inp log2ratio signal
#  -chip_vs_inp log2ratio signal
#  -chip subtracted log2ratio signal
#  -rzscore signal
#  -signed log p values
# For each of these, we generate the observed value plus a 95% CI.
# We assume that data are in hdf5 files.

if __name__ == "__main__":

    conf_file = sys.argv[1]
    conf_dict = toml.load(conf_file)
    conf_file_global = sys.argv[2]
    conf_dict_global = toml.load(conf_file_global)

    BINDIR = conf_dict_global["general"]["bindir"]
    PLOT = conf_dict_global["quant"]["diagnostic_plots"]
    ALPHA = conf_dict_global['quant']['alpha']
    MIN_PERC = conf_dict_global['quant']['min_percentile_chipsub_fit']
    SLOPE_THRESH = conf_dict_global['quant']['slope_increment_frac']

    # figure out some global parameters
    bs_opts = conf_dict_global['bootstrap']
    BS_NUM = bs_opts['bootstrap_samples']
    BS_DIR = bs_opts['bootstrap_direc']
    SAMP_TYPES = conf_dict["general"]["sample_types"]
    TYPE_COUNT = len(SAMP_TYPES)
    QNORM_TYPES = conf_dict["quant"]["qnorm_samples"]
    SPIKENORM_TYPES = conf_dict["quant"]["spikenorm_samples"]
    SPIKE_NAME = conf_dict_global['genome']['spike_in_name']

    # set up a lookup table to programatically associate sample types,
    # with their directories, file prefixes, and array indices later
    norm_lut = {
        'qnorm': {
            'dset' : conf_dict_global['norm']['qnorm_dset'],
            'type_lut' : {},
            'spikein_name' : SPIKE_NAME,
        },
        'spikenorm': {
            'dset' : conf_dict_global['norm']['spikenorm_dset'],
            'type_lut' : {},
            'spikein_name' : SPIKE_NAME,
        },
    }
    qnorm_idx = 0
    spikenorm_idx = 0
    for i,samp_type in enumerate(SAMP_TYPES):
        if samp_type in QNORM_TYPES:
            norm_lut['qnorm']['type_lut'][samp_type] = {
                'idx' : qnorm_idx,
                'direc' : conf_dict[samp_type]['directory'],
                'prefix' : conf_dict[samp_type]['sample_names']
            }
            qnorm_idx += 1
        if samp_type in SPIKENORM_TYPES:
            norm_lut['spikenorm']['type_lut'][samp_type] = {
                'idx' : spikenorm_idx,
                'direc' : conf_dict[samp_type]['directory'],
                'prefix' : conf_dict[samp_type]['sample_names']
            }
            spikenorm_idx += 1

    OUT_PREFIX = os.path.join(
        bs_opts['output_path'],
        conf_dict['general']['out_prefix']
    )

    #OUT_HDF_NAME = OUT_PREFIX + ".hdf5"

    # delete the old output hdf5 file if it already exists
    #if os.path.isfile(OUT_HDF_NAME):
    #    os.remove(OUT_HDF_NAME)

    write_outs = partial(
        qutils.write_outs,
        out_prefix = OUT_PREFIX,
        #out_hdf_name = OUT_HDF_NAME,
        #spike_name = SPIKE_NAME,
    )
    
    # list of samples to subtract chip from
    NUMER_LIST = conf_dict["quant"]["chipsub_numerators"]

    # make missing path if needed
    if not(os.path.isdir(bs_opts['output_path'])):
        os.mkdir(bs_opts['output_path'])

    # check for paired status. If user doesn't have paired
    #  in their configuration file, exit with an error instructing them
    #  to add the paired information to their config file.
    try:
        paired = conf_dict["quant"]["paired"]
    except KeyError:
        sys.exit("ERROR: You must include \"paired\" in your condition-level\
                  configuration file!! If you have unpiared data,\
                  int the [quant] block, you should create a new\
                  line that reads \"paired=false\", without the quotes.\
                  If you have paired data, write \"paired=true\", without\
                  the quotes in the [quant] block.")

    # NOTE: at some point it would be nice to be able to make this
    #   regexp term a configuration, but I can't figure out how
    #   to get yaml or toml to fully appreciate my brilliance....
    # Now we get our data array of shape (R,G,T) and modify type_lut
    #   in place to organize the data for future use. Specifically,
    #   type_lut['missing'] now exists, which contains a 2d array
    #   indicating which replicate_idx/type_idx are missing.
    regex_pat = re.compile(r'rep(\d+)')

    # norm_lut is modified in place here
    qutils.set_up_data_from_hdf2(
        norm_lut,
        conf_dict,
        BS_DIR,
        #pat = conf_dict["quant"]["rep_regexp"],
        pat = regex_pat,
    )

    #data_arr,missing_arr,ctg_lut,res = qutils.set_up_data_from_hdf(
    #    type_lut, # this dictionary is modified in place by this function
    #    conf_dict,
    #    BS_DIR,
    #    #pat = conf_dict["quant"]["rep_regexp"],
    #    pat = regex_pat,
    #    norm_dset_base = conf_dict_global['norm']['out_dset'],
    #)

    for norm_method,info in norm_lut.items():
        ctg_lut = info['ctg_lut']
        info['rev_ctg_lut'] = {
            ctg_info["idx"]: {
                "id": ctg_id, "length": ctg_info["length"]
            }
            for ctg_id,ctg_info in ctg_lut.items()
        }

    #hdf_utils.set_up_hdf_file(
    #    OUT_HDF_NAME,
    #    norm_lut['qnorm']['rev_ctg_lut'], #rev_ctg_lut,
    #    res,
    #    type_lut,
    #    paired,
    #)

    nprocs = conf_dict_global['quant']['quant_numproc']
    # If we have paired data then we probably don't need this many procs
    #   Reset nprocs to minimum of the current value or the number of
    #   chipsub calculations we're going to be doing later.
    #if paired:
    #    nprocs = int(np.min([nprocs,len(NUMER_LIST)*data_arr.shape[0]]))

    # If we have paired replicates, we impute every missing replicate
    #   within a given sample type here.
    # We impute using mean of the existing replicates for this sample
    #   type, and add noise. 
    #   The noise here came from bootstrapping in earlier steps.
    try:
        force = conf_dict['quant']['force_onesample_unpaired']
    except KeyError:
        force = False
    
    # Here we modify data_arr in place to supplement the missing 
    #   values.
    # Additionally, we modify missing_arr in place to switch the
    #   imputed R/T pairs to False if they were imputed here.
    for norm_method,info in norm_lut.items():
        
        spike_name = info['spikein_name']
        
        qutils.impute_missing_hdf(
            info['data_arr'],
            info['missing_arr'],
            info['type_lut'],
            BS_NUM,
            paired,
            force,
            #spike_name,
        )

    # At this point, if we don't have spike-in,
    #   we need to do median normalization to bring all
    #   of the various qnormed sample types' data into register.
    # The median_norm function modifies data in place to do just that
    # Default behavior is to set median for each replicate/sample type
    #   to 100.0
    # NOTE: if unpaired, there could still be replicate/sample type
    #   data for which all genome positions are zero. No worries,
    #   those indices will just be nan after median normalization.

    for norm_method,info in norm_lut.items():
        if norm_method == 'qnorm':
            qutils.median_norm(info['data_arr'])
    
    for norm_method,info in norm_lut.items():
        info['weights_arr'],info['jack_coefs'] = qutils.get_jackknife_repweights(
            info['data_arr'],
            info['missing_arr'],
            paired,
        )

    if paired:

        # If we have paired data, get the log2 ratios of each sample type
        #   vs input for each paired replicate. Save the result.
        # Then calculate chip-subtracted lot2 ratios on the desired data.
        # Proceed to jackknife sampling of log2 ratios and chipsub values.
        # Not providing a weights array causes the function to just compute
        #   the ratios within each replicate.
        for norm_method,info in norm_lut.items():
            info['log_rats'] = qutils.calc_lograt_vs_input(
                info['data_arr'],
                info['type_lut'],
            )
            write_outs(
                info['log_rats'],
                info['type_lut'],
                # here the sample type and replicate num will get substituted in
                #   within the write_outs2 function.
                info_str = '{}_vs_inp_lograt_rep{}',
                pat = regex_pat,
                #spike_name = info['spikein_name']
            )

            # Here we're calculating each jackknife replicate's log2_ratios
            #  for each sample type relative to input
            # jacked_log2_rats is shape (J,G,T). So, the first axis, J, contains
            #   the mean across sample replicates for each given jackknife
            #   repliate, j.
            info['jacked_log_rats'] = qutils.get_weighted_mean_within_jackknife_reps(
                info['log_rats'],
                info['weights_arr'],
            )

        (
            type_lut,
            log_rats,
            jacked_log_rats,
            weights_arr,
            jack_coefs,
            ctg_lut,
            rev_ctg_lut,
            res
        ) = qutils.gather_norm_data(
            norm_lut
        )

        # subtract trend in association between data of interest in NUMER_LIST
        #   and RNAP chip. This function does this in parallel, running nprocs
        #   processes at once.
        if NUMER_LIST: 

            chipsub_lut = qutils.get_chipsub_lut(
                type_lut,
                conf_dict["quant"]["chipsub_numerators"],
            )

            chipsub = qutils.do_chipsub(
                log_rats,
                type_lut,
                chipsub_lut,
                NUMER_LIST,
                plot_diagnostics = PLOT,
                chipsub_percentile = MIN_PERC,
                slope_threshold = SLOPE_THRESH,
                nproc = nprocs,
            )
            write_outs(
                chipsub,
                chipsub_lut,
                # here the sample type and replicate num will get substituted in
                #   within the write_outs2 function.
                info_str = '{}_chipsub_rep{}',
                pat = regex_pat,
            )

            # Calculate robust z-scores for each replicate's chipsub numbers.
            chipsub_rz = qutils.get_fn_over_axes(
                chipsub,
                iter_axis = [0,2],
                fn = qutils.calc_rzscores,
            )
            write_outs(
                chipsub_rz,
                chipsub_lut,
                info_str = "{}_rzchipsub_rep{}",
                pat = regex_pat,
            )
            # Calculate log10p from rz-scores for each replicate's chipsub.
            chipsub_log10p = qutils.get_fn_over_axes(
                chipsub_rz,
                iter_axis = [0,2],
                fn = qutils.calc_signed_log10p,
            )
            write_outs(
                chipsub_log10p,
                chipsub_lut,
                info_str = '{}_rzchipsublog10p_rep{}',
                pat = regex_pat,
            )
            # Now combine reps using jackknifing.
            # Get the sample type weights for appropriate numerators.
            orig_chipsub_idxs = [
                samp_info["orig_idx"]
                for _,samp_info in list(chipsub_lut.items())
            ]
            chipsub_weights = weights_arr[:,:,orig_chipsub_idxs]

            # Get the chipsub estimate for each jackknife replicate.        
            jacked_chipsub = qutils.get_weighted_mean_within_jackknife_reps(
                chipsub,
                chipsub_weights,
            )

    else:
        # Here we're calculating each jackknife replicate's log2_ratios
        #  for each normalization method and sample type relative to input
        # jacked_log_rats is shape (J,G,T).
        for norm_method,info in norm_lut.items():
            info['jacked_log_rats'] = qutils.calc_lograt_vs_input(
                info['data_arr'],
                info['type_lut'],
                info['weights_arr'],
            )

        type_lut,log_rats,jacked_log_rats = qutils.gather_norm_data(
            norm_lut
        )

        # jacked_chipsub is shape (J,G,N), where N is the number of numerators
        #   in numerator_list. chipsub_lut is a lookup table to associate sample
        #   types with appropriate indices in the N axis of jacked_chipsub.
        if NUMER_LIST: 

            chipsub_lut = qutils.get_chipsub_lut(
                type_lut,
                conf_dict["quant"]["chipsub_numerators"],
            )

            jacked_chipsub = qutils.do_chipsub(
                jacked_log_rats,
                type_lut,
                chipsub_lut,
                NUMER_LIST,
                plot_diagnostics = PLOT,
                chipsub_percentile = MIN_PERC,
                slope_threshold = SLOPE_THRESH,
                nproc = nprocs,
            )

    # Calculate robust z-scores for each jackknife replicate
    jacked_lograt_rz = qutils.get_fn_over_axes(
        jacked_log_rats,
        iter_axis = [0,2],
        fn = qutils.calc_rzscores,
    )
    # calculate jackknife-based mean, upper, lower conf limits
    lograt_rz_mean,lograt_rz_lo,lograt_rz_hi = qutils.calc_jackknife_cl(
        jacked_lograt_rz,
        jack_coefs,
        alpha = ALPHA,
    )

    # Calculate log10p robust z-scores for each jackknife replicate
    jacked_lograt_log10p = qutils.get_fn_over_axes(
        jacked_lograt_rz,
        iter_axis = [0,2],
        fn = qutils.calc_signed_log10p,
    )
    # calculate jackknife-based mean, uppoer, lower conf limits
    log10p_lr_mean,log10p_lr_lo,log10p_lr_hi = qutils.calc_jackknife_cl(
        jacked_lograt_log10p,
        jack_coefs,
        alpha = ALPHA,
    )


    # calculate mean estimate and conf lims accross jackknife replicates
    log_mean,log_lo,log_hi = qutils.calc_jackknife_cl(
        jacked_log_rats,
        jack_coefs,
        alpha = ALPHA,
    )

    # calculate jackknife-based mean, upper, lower conf limits
    if NUMER_LIST: 
        chipsub_mean,chipsub_lo,chipsub_hi = qutils.calc_jackknife_cl(
            jacked_chipsub,
            jack_coefs,
            alpha = ALPHA,
        )

        # Calculate robust z-scores for each jackknife replicate
        jacked_chipsub_rz = qutils.get_fn_over_axes(
            jacked_chipsub,
            iter_axis = [0,2],
            fn = qutils.calc_rzscores,
        )
        # calculate jackknife-based mean, upper, lower conf limits
        chipsub_rz_mean,chipsub_rz_lo,chipsub_rz_hi = qutils.calc_jackknife_cl(
            jacked_chipsub_rz,
            jack_coefs,
            alpha = ALPHA,
        )

        # Calculate log10p robust z-scores for each jackknife replicate
        jacked_chipsub_log10p = qutils.get_fn_over_axes(
            jacked_chipsub_rz,
            iter_axis = [0,2],
            fn = qutils.calc_signed_log10p,
        )
        # calculate jackknife-based mean, uppoer, lower conf limits
        log10p_mean,log10p_lo,log10p_hi = qutils.calc_jackknife_cl(
            jacked_chipsub_log10p,
            jack_coefs,
            alpha = ALPHA,
        )

    write_outs(
        log_mean,
        type_lut,
        info_str = '{}_vs_inp_lograt_mean',
    )
    write_outs(
        log_lo,
        type_lut,
        info_str = '{{}}_vs_inp_lograt_{}cllo'.format(int(ALPHA*100)),
    )
    write_outs(
        log_hi,
        type_lut,
        info_str = '{{}}_vs_inp_lograt_{}clhi'.format(int(ALPHA*100)),
    )
    write_outs(
        lograt_rz_mean,
        type_lut,
        info_str = '{}_vs_inp_rzlograt_mean',
    )
    write_outs(
        lograt_rz_lo,
        type_lut,
        info_str = '{{}}_vs_inp_rzlograt_{}cllo'.format(int(ALPHA*100)),
    )
    write_outs(
        lograt_rz_hi,
        type_lut,
        info_str = '{{}}_vs_inp_rzlograt_{}clhi'.format(int(ALPHA*100)),
    )
    write_outs(
        log10p_lr_mean,
        type_lut,
        info_str = '{}_vs_inp_rzlogratlog10p_mean',
    )
    write_outs(
        log10p_lr_lo,
        type_lut,
        info_str = '{{}}_vs_inp_rzlogratlog10p_{}cllo'.format(int(ALPHA*100)),
    )
    write_outs(
        log10p_lr_hi,
        type_lut,
        info_str = '{{}}_vs_inp_rzlogratlog10p_{}clhi'.format(int(ALPHA*100)),
    )

    if NUMER_LIST: 
        write_outs(
            chipsub_mean,
            chipsub_lut,
            info_str = '{}_chipsub_mean',
        )
        write_outs(
            chipsub_lo,
            chipsub_lut,
            info_str = '{{}}_chipsub_{}cllo'.format(int(ALPHA*100)),
        )
        write_outs(
            chipsub_hi,
            chipsub_lut,
            info_str = '{{}}_chipsub_{}clhi'.format(int(ALPHA*100)),
        )
        write_outs(
            chipsub_rz_mean,
            chipsub_lut,
            info_str = '{}_rzchipsub_mean'
        )
        write_outs(
            chipsub_rz_lo,
            chipsub_lut,
            info_str = "{{}}_rzchipsub_{}cllo".format(int(ALPHA*100)),
        )
        write_outs(
            chipsub_rz_hi,
            chipsub_lut,
            info_str = "{{}}_rzchipsub_{}clhi".format(int(ALPHA*100)),
        )
        write_outs(
            log10p_mean,
            chipsub_lut,
            info_str = '{}_rzchipsublog10p_mean'
        )
        write_outs(
            log10p_lo,
            chipsub_lut,
            info_str = "{{}}_rzchipsublog10p_{}cllo".format(int(ALPHA*100)),
        )
        write_outs(
            log10p_hi,
            chipsub_lut,
            info_str = "{{}}_rzchipsublog10p_{}clhi".format(int(ALPHA*100)),
        )

