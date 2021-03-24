#!/usr/bin/python

# calculate peak calls for all of my ipod samples of interest
# we read all of the instances to look at from a table of name/config file pairs

import argparse
import os
import os.path
import toml
import subprocess


# the following command takes three arguments: an input .gr file, an output .gr file, and a threshold value for peak calls
PEAK_CALL_SCRIPT="python /data/petefred/st_lab_work/ipod_hr/bin/fire/call_peaks_v7.py %s %s %i %f"

# this one just need the peaks .gr file and output prefix
OVERLAP_SCRIPT= "python /home/petefred/src/pytools/analyze_peaks.py %s /data/petefred/st_lab_work/e_coli_data/regulondb_20180516/BindingSites_knownsites_flags.gr > %s_tf_overlaps.txt"

# set up to parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument('--basedir',help="root directory to find all files of interest (default: current directory)", default=os.getcwd())
parser.add_argument('--outdir',help="root directory for output files (default: current directory)", default=os.getcwd())
parser.add_argument('--windowsize',help="Window size for rolling average (default 15)", type=int, default=15)
parser.add_argument('cond_list', help='File containing a list of condition directories and config files', metavar='condition_list')

args=parser.parse_args()

# now go through the conditions of interest and run the analysis
# we actually call the peaks, and then compare them to tfbs lists

conf_str=open(args.cond_list)
for line in conf_str:
    print line
    dirname,conffile = line.rstrip().split()
    cond_conf = toml.load(os.path.join(args.basedir, dirname, conffile))
    gr_file = os.path.join(args.basedir, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + "_v6rz_chipsub.gr")
    for cutoff in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
        run_cmd = PEAK_CALL_SCRIPT % ( gr_file, os.path.join(args.outdir, dirname + "_cutoff%f_peaks.gr" % cutoff), args.windowsize, cutoff)
        subprocess.call(run_cmd, shell=True)
        analyze_cmd = OVERLAP_SCRIPT % (os.path.join(args.outdir, dirname + "_cutoff%f_peaks.gr" % cutoff), os.path.join(args.outdir, dirname + "_cutoff%f_peaks.gr" % cutoff))
        subprocess.call(analyze_cmd, shell=True)

    gr_file = os.path.join(args.basedir, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + "_v6rzlog10p_chipsub.gr")
    for cutoff in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0]:
        run_cmd = PEAK_CALL_SCRIPT % ( gr_file, os.path.join(args.outdir, dirname + "_log10p_cutoff%f_peaks.gr" % cutoff), args.windowsize, cutoff)
        subprocess.call(run_cmd, shell=True)
        analyze_cmd = OVERLAP_SCRIPT % (os.path.join(args.outdir, dirname + "_log10p_cutoff%f_peaks.gr" % cutoff), os.path.join(args.outdir, dirname + "_log10p_cutoff%f_peaks.gr" % cutoff))
        subprocess.call(analyze_cmd, shell=True)
