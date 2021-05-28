#!/usr/bin/python

# script to run fastqc on each input file in my sequencing manifest
# these should be inspected to make sure everything looks ok, likely before doing ANYTHING else
# by default we run this on the clipped/trimmed reads that are ready to be aligned

import subprocess
import os
import toml
import sys

conf_file = sys.argv[1]
conf_dict = toml.load(conf_file)
conf_file_global = sys.argv[2]
conf_dict_global = toml.load(conf_file_global)
samp_type = sys.argv[3]

# useful constants
READDIR = conf_dict_global["processing"]["processed_direc"] # directory containing the trimmed reads and fastqc output
THREADS = conf_dict_global["qc"]["fastqc_threads"]
QCDIR = conf_dict_global["qc"]["qc_direc"]

F_SUFFIX = conf_dict_global["processing"]["f_paired_read_file_suffix"]
R_SUFFIX = conf_dict_global["processing"]["r_paired_read_file_suffix"]

def run_fastqc_fqfiles(infile_f, infile_r, n_threads):
    outdir = os.path.join(QCDIR, 'fastqc_after')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    cmdline = "fastqc {} {} --noextract -t {} -o {}".format(
        infile_f, infile_r, n_threads, outdir
    )
    subprocess.call(cmdline,shell=True)

samp_prefixes = conf_dict[samp_type][sample_names]

for inprefix in samp_prefixes:
    infile_fwd = os.path.join(READDIR, inprefix + F_SUFFIX)
    infile_rev = os.path.join(READDIR, inprefix + R_SUFFIX)
    run_fastqc_fqfiles(infile_fwd, infile_rev, THREADS)

