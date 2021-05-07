#!/usr/bin/python

# calculate peak calls for all ipod samples of interest
# we read all of the instances to look at from a table of name/config file pairs

import argparse
import os
import sys
import toml
import subprocess
import glob
import re
import numpy as np

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../src/utils")
sys.path.insert(0, utils_path)

import hdf_utils
import anno_tools as anno

# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    'main_conf',
    help="Configuration file defining work to be done.",
)
parser.add_argument(
    '--skipsteps',
    help="comma-separated list of steps to skip. Can be any of (peaks,epods)."
)
args = parser.parse_args()

if args.skipsteps is None:
    skipsteps = set()
else:
    skipsteps = set(args.skipsteps.split(','))

steps = ['peaks','epods','idr']
for step in skipsteps:
    if step not in steps:
        sys.exit("\nERROR: {} is not a step. Allowed steps are peaks, epods, idr.\n".format(step))

# parse the top level config file to get some needed information
conf_file = args.main_conf
conf_dict_global = toml.load(conf_file)

BASEDIR = conf_dict_global["general"]["basedir"]
BINDIR = conf_dict_global["general"]["bindir"]
WINSIZE = conf_dict_global["peaks"]["windowsize"]
SAMP_FNAME = os.path.join(
    BASEDIR,
    conf_dict_global["general"]["condition_list"],
)

RESOLUTION = conf_dict_global["genome"]["resolution"]
SEQ_DB = conf_dict_global["genome"]["genome_base"]

# the following command takes three arguments: an input .gr file, an output .gr file, and a threshold value for peak calls
PEAK_CALL_SCRIPT = "python {}/peakcalling/call_peaks.py\
                        --in_file {{}}\
                        --sample_type {{}}\
                        --out_file {{}}\
                        --window_size {}\
                        --threshold {{}}".format(BINDIR,WINSIZE)

PEAK_IDR_SCRIPT = "idr --samples {} {}\
                       --plot --log-output-file {}.log --verbose\
                       --output-file {}"

# this one just need the peaks .gr file and output prefix
OVERLAP_SCRIPT = "python {}/peakcalling/analyze_peaks.py {{}}\
                  /data/petefred/st_lab_work/e_coli_data/regulondb_20180516/BindingSites_knownsites_flags.gr > \
                 {{}}_tf_overlaps.txt".format(
    BINDIR,
)

EPOD_CALL_SCRIPT = "python {}/epodcalling/call_epods.py\
                        --in_file {{}}\
                        --sample_type {{}}\
                        --out_file {{}}".format(BINDIR)

## get contig lengths using hdf_utils.make_ctg_lut_from_bowtie
## then make arrays for each contig to store peak loci passing
## IDR threshold
ctg_lut = hdf_utils.make_ctg_lut_from_bowtie(SEQ_DB)
ctg_array_dict = {}
for ctg_idx,ctg_info in ctg_lut.items():
    ctg_len = ctg_info["length"]
    # now we have a dictionary with ctg id as keys, zeros array as vals
    ctg_array_dict[ctg_info["id"]] = {}
    ctg_array_dict[ctg_info["id"]]["loci"] = np.arange(0, ctg_len, RESOLUTION)

# now go through the conditions of interest and run the analysis
# we actually call the peaks, and then compare them to tfbs lists

samp_file = open(SAMP_FNAME)
for line in samp_file:

    dirname,samp_conf = line.rstrip().split()
    dir_path = os.path.join(BASEDIR, dirname)
    os.chdir(dir_path)
    conf_dict = toml.load(os.path.join(dir_path, samp_conf))
    chipsub_samps = conf_dict["quant"]["chipsub_numerators"]
    no_chipsub_samps = conf_dict["quant"]["no_chipsub"]
    out_path = os.path.join(dir_path, conf_dict["general"]["output_path"])
    paired = conf_dict["quant"]["paired"]

    all_samps = []
    all_samps.extend(chipsub_samps)
    all_samps.extend(no_chipsub_samps)

    cutoff_dict = {
        'rz': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        'log10p': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                   10.0, 12.0, 15.0, 20.0, 30.0, 50.0],
    }

    for score_type in ['rz','log10p']:

        # loop over all samples
        for samp in all_samps:

            # if the sample was in the chipsub category, its dset name
            #   looks something like this
            if samp in chipsub_samps:
                if score_type == 'rz':
                    fname = "{}_{}_rzchipsub_{{}}.bedgraph".format(
                        dirname, samp.upper()
                    )
                elif score_type == 'log10p':
                    fname = "{}_{}_rzchipsublog10p_{{}}.bedgraph".format(
                        dirname, samp.upper()
                    )

            # if the sample did NOT have chipsub performed, its dset name
            #   looks something like this.
            else:
                if score_type == 'rz':
                    fname = "{}_{}_vs_inp_rzlograt_{{}}.bedgraph".format(
                        dirname, samp.upper()
                    )
                elif score_type == 'log10p':
                    fname = "{}_{}_vs_inp_rzlogratlog10p_{{}}.bedgraph".format(
                        dirname, samp.upper()
                    )

            # If these data were not from paired samples of inp/chip/ipod,
            #   then just use the mean result for peak calling
            if not paired:
                fname_list = [ fname.format("mean") ]
            # If the data were from paired samples of inp/chip/ipod,
            #   then get each replicate's dataset name.
            else:
                fname_search = fname.format("rep*")
                fname_list = glob.glob(os.path.join(out_path, fname_search))

            # do peak calling
            if not 'peaks' in skipsteps:
                # loop over multiple score cutoffs.
                for cutoff in cutoff_dict[score_type]:

                    out_files = []

                    # loop over files. Just one if it's not paired data.
                    for fname in fname_list:

                        base_name = os.path.basename(fname)
                        base_name_prefix = os.path.splitext(base_name)[0]

                        out_np_path = os.path.join(
                            out_path,
                            "{}_cutoff_{}_peaks.narrowpeak".format(
                                base_name_prefix,
                                cutoff,
                            ),
                        )

                        run_cmd = PEAK_CALL_SCRIPT.format(
                            fname,
                            samp,
                            out_np_path,
                            cutoff,
                        )
                        subprocess.call(run_cmd, shell=True)
                        
                        out_files.append(out_np_path)

                    if not 'idr' in skipsteps:
                        if paired:
                            # go over replicates' peaks and do pair-wise IDR calculation
                            #   for each pair-wise grouping of replicates
                            # Save narrowpeak output for each IDR calculation
                            rep_count = len(out_files)
                            n_idrs = (rep_count**2 - rep_count) / 2
                            for ctg_idx,ctg_info in ctg_lut.items():
                                ctg_len = ctg_info["length"]
                                ctg_array_dict[ctg_info["id"]]["num_passed_array"] = np.zeros(ctg_len)

                            rep_idxs = np.asarray([i for i in range(rep_count)])
                            idr_outfiles = []
                            
                            for idx_a in rep_idxs:
                                for idx_b in rep_idxs[rep_idxs > idx_a]:

                                    fname_a = out_files[idx_a]
                                    pref_a = os.path.splitext(
                                        os.path.basename(fname_a)
                                    )[0]
                                    fname_b = out_files[idx_b]
                                    pref_b = os.path.splitext(
                                        os.path.basename(fname_b)
                                    )[0]
                                    idr_out_pref = "{}_vs_{}_idr".format(
                                        pref_a,
                                        pref_b,
                                    )
                                    idr_out_pref = os.path.join(out_path, idr_out_pref)
                                    idr_outfile = idr_out_pref + ".narrowpeak"
                                    idr_outfiles.append(idr_outfile)
                                    
                                    print("Calculating IDR for each peak in {} and {}.".format(fname_a, fname_b))
                                    idr_cmd = PEAK_IDR_SCRIPT.format(
                                        fname_a,
                                        fname_b,
                                        idr_out_pref,
                                        idr_outfile,
                                    )
                                    subprocess.call(idr_cmd, shell=True)

                            # iterate over idr comparisons and add 1 to each
                            #   position passing IDR < 0.05
                            for idx,idr_file in enumerate(idr_outfiles):
                                idr_results = anno.NarrowPeakData()
                                idr_results.parse_narrowpeak_file(idr_file)
                                
                                # each peak in the idr file has a global idr we
                                #   filter by, then if passing that filter,
                                #   add 1 to the genomic region encompassed
                                #   by this peak.
                                for peak in idr_results:
                                    if peak.global_idr < 0.05:
                                        ctg_array_dict[peak.chrom_name][
                                            "num_passed_array"
                                        ][
                                            peak.start:peak.end
                                        ] += 1
                            frac_passed_dict = {}
                            
                            # caluclate fraction of comparisons that are 
                            #   reproducible for this threshold
                            for ctg_id,ctg_info in ctg_array_dict.items():
                                frac_passed_arr = (
                                    ctg_info["num_passed_array"] / len(idr_outfiles)
                                )
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################

            # do epod calling
            if not 'epods' in skipsteps:
                
                out_np_path = os.path.join(
                    out_path,
                    "{}_{}_epods.narrowpeak".format(prefix,dset),
                )

                epod_cmd = EPOD_CALL_SCRIPT.format(
                    in_bg_path,
                    samp,
                    dset_name,
                    out_np_path,
                )
                subprocess.call(epod_cmd, shell=True)

    #        analyze_cmd = OVERLAP_SCRIPT.format(
    #            os.path.join(
    #                args.outdir,
    #                dirname + "_rz_cutoff_{}_peaks.gr".format(cutoff),
    #            ),
    #            os.path.join(
    #                args.outdir,
    #                dirname + "_rz_cutoff_{}_peaks.gr".format(cutoff),
    #            ),
    #        )
    #        subprocess.call(analyze_cmd, shell=True)
    #
    #    gr_file = os.path.join(
    #        args.basedir,
    #        dirname,
    #        conf_dict["general"]["output_path"],
    #        conf_dict["general"]["out_prefix"] + "_v6rzlog10p_chipsub.gr",
    #    )
    #
    #    for cutoff in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0]:
    #
    #        run_cmd = PEAK_CALL_SCRIPT.format(
    #            gr_file,
    #            os.path.join(
    #                args.outdir,
    #                dirname + "_log10p_cutoff_{}_peaks.gr".format(cutoff),
    #            ),
    #            cutoff,
    #        )
    #        subprocess.call(run_cmd, shell=True)
    #
    #        analyze_cmd = OVERLAP_SCRIPT.format(
    #            os.path.join(
    #                args.outdir,
    #                dirname + "_log10p_cutoff_{}_peaks.gr".format(cutoff),
    #            ),
    #            #NOTE: I think we want to get rid of the .gr suffix below.
    #            os.path.join(
    #                args.outdir,
    #                dirname + "_log10p_cutoff_{}_peaks.gr".format(cutoff),
    #            ),
    #        )
    #        subprocess.call(analyze_cmd, shell=True)

