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

steps = ['peaks','epods']
for step in skipsteps:
    if step not in steps:
        sys.exit("\nERROR: {} is not a step. Allowed steps are peaks, epods.\n".format(step))

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

# the following command takes three arguments: an input .gr file, an output .gr file, and a threshold value for peak calls
PEAK_CALL_SCRIPT = "python {}/peakcalling/call_peaks.py\
                        --in_file {{}}\
                        --sample_type {{}}\
                        --out_file {{}}\
                        --window_size {}\
                        --threshold {{}}".format(BINDIR,WINSIZE)

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

                    ####################################################
                    ####################################################
                    ######## here I'll set up an array for results of peak calling
                    ######## it can then be used for idr ###############
                    ####################################################
                    ####################################################



                    # loop over files. Just one if it's not paired data.
                    for fname in fname_list:

                        base_name = os.path.basename(fname)
                        base_name_prefix = os.path.splitext(base_name)[0]

                        out_bg_path = os.path.join(
                            out_path,
                            "{}_cutoff_{}_peaks.narrowpeak".format(
                                base_name_prefix,
                                cutoff,
                            ),
                        )

                        run_cmd = PEAK_CALL_SCRIPT.format(
                            fname,
                            samp,
                            out_bg_path,
                            cutoff,
                        )
                        subprocess.call(run_cmd, shell=True)

                        

            # do epod calling
            if not 'epods' in skipsteps:
                
                out_bg_path = os.path.join(
                    out_path,
                    "{}_{}_epods.bedgraph".format(prefix,dset),
                )

                epod_cmd = EPOD_CALL_SCRIPT.format(
                    in_bg_path,
                    samp,
                    dset_name,
                    out_bg_path,
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

