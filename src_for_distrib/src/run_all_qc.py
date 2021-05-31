#!/usr/bin/python

# Run quality control for all samples in a given config file.
# This is intended to be run from the top level directory for those samples.
# The config file should indicate in the [general]->sample_types option
#   what types of samples (inp, chip, etc) are present; there
#   should already be a directory for each such sample type.
# In addition, each of those directories should have a properly filled out
#   read_manifest.txt, and must have raw/, aligned/, and bootstrap/ directories
#   are already set up.

import sys
import subprocess
import toml
import os

# at this point, if we've run run_all_driver.py, we've 
# now moved into a data directory, such as wt_0h
conf_file = sys.argv[1] # Abs path to condition-level conf file
conf_dict = toml.load(conf_file)
conf_file_global = sys.argv[2] # Abs path to global configuration here
conf_dict_global = toml.load(conf_file_global)

# get location of scripts
SRC_DIR = conf_dict_global["general"]["bindir"]
BSDIR = conf_dict_global["bootstrap"]["bootstrap_direc"]

QC_CMD = "python {}/run_fastqc_onedir.py {} {} {{}}".format(
    SRC_DIR, conf_file, conf_file_global
)
CHIP_CMD = "python {}/run_chipqc_onedir.py {} {{}}".format(
    SRC_DIR, conf_file_global
)
BASE_DIR = os.getcwd()

def run_qc_driver(hdf, sampname, sampdir, n_errors):
    # helper function to run qc on a specified sample type
    print("Working on {} samples...".format(sampname))
    try:
        os.chdir(sampdir)
        subprocess.check_call(QC_CMD.format(sampname), shell=True)
        subprocess.check_call(CHIP_CMD.format(hdf), shell=True)
        os.chdir(BASE_DIR)
    except:
        n_errors += 1
        print("Warning: Encountered an error processing {} samples".format(sampname))

    return n_errors

n_errors = 0

print("Now running fastqc on all preprocessed data...")

# for each of ["inp","chip","ipod",...], run qc
for sample_type in conf_dict["general"]["sample_types"]:
    
    sample_prefixes = conf_dict[sample_type]["sample_names"]
    sample_dir = conf_dict[sample_type]["directory"]

    hdf_names = [
        os.path.join( BSDIR, pref ) + ".hdf5"
        for pref in sample_prefixes
    ]

    for hdf_name in hdf_names:
        n_errors = run_qc_driver(
            hdf_name,
            sample_type,
            conf_dict[sample_type]["directory"],
            n_errors
        )

print("Finished with all QC. Encountered {} errors".format(n_errors))
