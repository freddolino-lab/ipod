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
import multiprocessing

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../src/utils")
sys.path.insert(0, utils_path)

import hdf_utils
import anno_tools as anno
import peak_utils as pu

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
parser.add_argument(
    '--invert_scores',
    help="Settin this option will call extended regions of depleted protein occupancy",
    action="store_true",
)
args = parser.parse_args()

if args.skipsteps is None:
    skipsteps = set()
else:
    skipsteps = set(args.skipsteps.split(','))

steps = ['peaks','epods']
for step in skipsteps:
    if step not in steps:
        sys.exit("\nERROR: {} is not a step. Allowed steps are peaks and epods.\n".format(step))

# parse the top level config file to get some needed information
conf_file = args.main_conf
conf_dict_global = toml.load(conf_file)

BASEDIR = conf_dict_global["general"]["basedir"]
BINDIR = conf_dict_global["general"]["bindir"]
RESOLUTION = conf_dict_global["genome"]["resolution"]
WINSIZE = int(conf_dict_global["peaks"]["windowsize_bp"] / RESOLUTION)
SAMP_FNAME = os.path.join(
    BASEDIR,
    conf_dict_global["general"]["condition_list"],
)
SEQ_DB = conf_dict_global["genome"]["genome_base"]
PEAK_PROCS = conf_dict_global["peaks"]["nproc"]
EPOD_PROCS = conf_dict_global["epods"]["nproc"]

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

## this one just need the peaks .gr file and output prefix
#OVERLAP_SCRIPT = "python {}/peakcalling/analyze_peaks.py {{}}\
#                  /data/petefred/st_lab_work/e_coli_data/regulondb_20180516/BindingSites_knownsites_flags.gr > \
#                 {{}}_tf_overlaps.txt".format(
#    BINDIR,
#)

EPOD_CALL_SCRIPT = "python {}/epodcalling/call_epods.py\
                        --main_conf {}\
                        --in_file {{}}\
                        --out_prefix {{}}\
                        --resolution {}".format(
    BINDIR,
    os.path.abspath(args.main_conf),
    RESOLUTION,
)

if args.invert_scores:
    EPOD_CALL_SCRIPT += " --invert_scores"

# read in toml file containing info on singularity versions if we're running this
#   from within a singularity container
if "IPOD_VER" in os.environ:
    VERSION = os.environ["IPOD_VER"]
    ver_filepath = os.path.join(BASEDIR, "singularity_version_info.toml")
    if os.path.isfile(ver_filepath):
        ver_info = toml.load(ver_filepath)
    else:
        ver_info = {
            "preprocessing": None,
            "alignment": None,
            "bootstrapping": None,
            "qc": None,
            "qnorm": None,
            "quant": None,
            "peak_calls": None,
            "epod_calls": None,
        }


def call_peaks(in_fname, out_path, cutoff, samp_name):
    '''Utility function used to wrap creation of peak calling subprocess
    into a function, thus allowing simple addition to a multiprocessing
    pool.

    Args:
    -----
    in_fname : str
        Input file name
    out_path : str
        Path to which results should be written
    cutoff : float
        Value above which regions will be called peaks
    samp_name : str
        Name of this sample, i.e., IPOD, NC, etc.

    Returns:
    --------
    Path to the narrowpeak file generated by peak calling.
    '''
    base_name = os.path.basename(in_fname)

    base_name_prefix = os.path.splitext(base_name)[0]

    out_np_path = os.path.join(
        out_path,
        "{}_cutoff_{}_peaks.narrowpeak".format(
            base_name_prefix,
            cutoff,
        ),
    )

    run_cmd = PEAK_CALL_SCRIPT.format(
        in_fname,
        samp_name,
        out_np_path,
        cutoff,
    )
    subprocess.call(run_cmd, shell=True)

    return out_np_path

def call_epods(in_fname, out_path):
    '''Utility function used to wrap creation of epod calling subprocess
    into a function, thus allowing simple addition to a multiprocessing
    pool.

    Args:
    -----
    in_fname : str
        Input file name
    out_path : str
        Path to which results should be written

    Returns:
    --------
    Path to the narrowpeak file generated by epod calling.
    '''

    base_name = os.path.basename(in_fname)
    base_name_prefix = os.path.splitext(base_name)[0]

    out_prefix = os.path.join(
        out_path,
        base_name_prefix,
    )
    epod_outfile = out_prefix + "_epods.narrowpeak"
    strict_epod_outfile = out_prefix + "_epods_strict.narrowpeak"

    epod_cmd = EPOD_CALL_SCRIPT.format(
        in_fname,
        out_prefix,
    )
    subprocess.call(epod_cmd, shell=True)

    return (epod_outfile, strict_epod_outfile)

def generate_fname(samp, chipsub_samps, score_type, out_prefix):

    # if the sample was in the chipsub category, its dset name
    #   looks something like this
    if samp in chipsub_samps:
        if score_type == 'rz':
            fname = "{}_{}_rzchipsub_{{}}.bedgraph".format(
                out_prefix, samp.upper()
            )
        elif score_type == 'log10p':
            fname = "{}_{}_rzchipsublog10p_{{}}.bedgraph".format(
                out_prefix, samp.upper()
            )

    # if the sample did NOT have chipsub performed, its dset name
    #   looks something like this.
    else:
        if score_type == 'rz':
            fname = "{}_{}_vs_inp_rzlograt_{{}}.bedgraph".format(
                out_prefix, samp.upper()
            )
        elif score_type == 'log10p':
            fname = "{}_{}_vs_inp_rzlogratlog10p_{{}}.bedgraph".format(
                out_prefix, samp.upper()
            )

    return fname

def calc_idr(paired, out_files, ctg_lut, out_path, fname, mean_fname, in_path, idr_thresh, cutoff=None):

    if paired:
        # go over replicates' peaks and do pair-wise IDR calculation
        #   for each pair-wise grouping of replicates
        # Save narrowpeak output for each IDR calculation
        rep_count = len(out_files)
        n_idrs = (rep_count**2 - rep_count) / 2
        for ctg_idx,ctg_info in ctg_lut.items():
            ctg_len = ctg_info["length"]
            ctg_array_dict[ctg_info["id"]]["num_passed_array"] = np.zeros(int(ctg_len/RESOLUTION))

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
                idr_out_pref = os.path.join(
                    out_path,
                    idr_out_pref
                )
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

        base_name = os.path.basename(fname)
        pu.compile_idr_results(
            idr_outfiles,
            ctg_array_dict,
            RESOLUTION,
            base_name,
            mean_fname,
            in_path,
            out_path,
            idr_thresh,
            cutoff = cutoff,
        )

def process_sample(line, conf_dict_global):

    dirname,samp_conf = line.rstrip().split()
    dir_path = os.path.join(BASEDIR, dirname)
    os.chdir(dir_path)
    conf_dict = toml.load(os.path.join(dir_path, samp_conf))
    out_file_prefix = conf_dict["general"]["out_prefix"]
    chipsub_samps = conf_dict["quant"]["chipsub_numerators"]
    no_chipsub_samps = conf_dict["quant"]["no_chipsub"]

    in_path = os.path.join(dir_path, conf_dict_global["bootstrap"]["output_path"])
    peak_out_path = os.path.join(dir_path, conf_dict_global["peaks"]["output_path"])
    epod_out_path = os.path.join(dir_path, conf_dict_global["epods"]["output_path"])

    if not os.path.isdir(peak_out_path):
        os.mkdir(peak_out_path)
    if not os.path.isdir(epod_out_path):
        os.mkdir(epod_out_path)

    paired = conf_dict["quant"]["paired"]

    all_samps = []
    all_samps.extend(chipsub_samps)
    all_samps.extend(no_chipsub_samps)

    cutoff_dict = {
        'rz': conf_dict_global["peaks"]["rz_thresholds"],
        'log10p': conf_dict_global["peaks"]["log10p_thresholds"],
    }

    idr_threshold = conf_dict_global["idr"]["threshold"]

    for score_type in ['rz','log10p']:

        # loop over all samples
        for samp in all_samps:

            fname_base = generate_fname(
                samp,
                chipsub_samps,
                score_type,
                out_file_prefix,
            )

            fname = os.path.join(in_path, fname_base)

            # If these data were not from paired samples of inp/chip/ipod,
            #   then just use the mean result for peak calling
            if not paired:
                fname_list = [ fname.format("mean") ]
                mean_fname = fname_list[0]
            # If the data were from paired samples of inp/chip/ipod,
            #   then get each replicate's dataset name.
            else:
                fname_search = fname.format("rep*")
                fname_list = glob.glob(os.path.join(in_path, fname_search))
                mean_fname = fname.format("mean")

            # do peak calling
            if not 'peaks' in skipsteps:
                # loop over multiple score cutoffs.
                for cutoff in cutoff_dict[score_type]:

                    # loop over files. Just one if it's not paired data.
                    out_files = []
                    for peak_fname in fname_list:

                        out_files.append(
                            call_peaks(
                                peak_fname,
                                peak_out_path,
                                cutoff,
                                samp,
                            )
                        )

                    if len(fname_list) < 2:
                        print("============================")
                        print("Only one replicate found in sample {}, skipping IDR calculation for peaks.".format(samp))
                        print("----------------------------")
                        continue

                    calc_idr(
                        paired,
                        out_files,
                        ctg_lut,
                        peak_out_path,
                        fname,
                        mean_fname,
                        in_path,
                        idr_threshold,
                        cutoff,
                    )

            # do epod calling
            if not 'epods' in skipsteps:

                if score_type == 'log10p':
                    continue

                strict_epod_outfiles = []
                loose_epod_outfiles = []

                for fname in fname_list:

                    these_outfiles = call_epods(
                        fname,
                        epod_out_path,
                    )
                    loose_epod_outfiles.append(these_outfiles[0])
                    strict_epod_outfiles.append(these_outfiles[1])

                if len(fname_list) < 2:
                    print("============================")
                    print("Only one replicate found in sample {}, skipping IDR calculation for EPODs.".format(samp))
                    print("----------------------------")
                    continue

                calc_idr(
                    paired,
                    loose_epod_outfiles,
                    ctg_lut,
                    epod_out_path,
                    fname,
                    mean_fname,
                    in_path,
                    idr_threshold,
                )
                calc_idr(
                    paired,
                    strict_epod_outfiles,
                    ctg_lut,
                    epod_out_path,
                    fname,
                    mean_fname,
                    in_path,
                    idr_threshold,
                )

    return fname


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

# use multiprocessing to do all of this in parallel
pool = multiprocessing.Pool(EPOD_PROCS)
all_res = []

samp_file = open(SAMP_FNAME)
for line in samp_file:

    all_res.append(
        pool.apply_async(
            process_sample,
            [
                line,
                conf_dict_global,
            ]
        )
    )

pool.close()
pool.join()

n_err = 0
for res in all_res:
    if not res.successful():
        print("\n==============================")
        print("Encountered error processing {}.".format(res.get()))
        print("------------------------------\n")
        n_err += 1

print("Finished running peak and epod calling jobs. Encountered {} errors.".format(n_err))

if "IPOD_VER" in os.environ:
    if n_err == 0:
        ver_info["peak_calls"] = VERSION
        ver_info["epod_calls"] = VERSION
        with open(ver_filepath, "w") as f:
            toml.dump(ver_info, f)

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

