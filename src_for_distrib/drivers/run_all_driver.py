#!/usr/bin/python

# actually run all of the needed commands in this directory

import subprocess
import sys
import os 
import toml
import argparse
import numpy as np


# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    'main_conf',
    help="main configuration file defining work to be done"
)
parser.add_argument(
    '--skipsteps',
    help="comma-separated list of steps to skip\
         (can be any of preprocess,align,bootstrap,qc,qnorm,quant)",
    default=None
)

args = parser.parse_args()
conf_file = args.main_conf

if args.skipsteps is None:
    skipsteps = set()
else:
    skipsteps = set(args.skipsteps.split(","))

steps = ['preprocess', 'align','bootstrap','qc','qnorm','spikenorm','quant']
for step in skipsteps:
    if step not in steps:
        sys.exit("\nERROR: {} is not an allowable step to skip. Allowed steps are preprocess, align, bootstrap, qc, qnorm, spikenorm, quant.\n".format(step))

# parse the top level config file to get some needed information
conf_dict_global = toml.load(conf_file)

BASEDIR = conf_dict_global["general"]["basedir"]
BINDIR = conf_dict_global["general"]["bindir"]
SPIKE = False
SPIKE_CHR = conf_dict_global["genome"]["spike_in_name"]
if SPIKE_CHR != "None":
    SPIKE = True
    # If we are doing spike-in, we'll skip quantile normalization
    #  so add qnorm to skipsteps
    skipsteps.add('qnorm')

# if we've run this driver from within our singularity container,
#   then the IPOD_VER environment veriable with exist.
# In that case, read in the toml file (if it exists) denoting which versions of
#   the container were used for which steps in the past. If the file doesn't 
#   exist, create a dictionary populated with correct information, and save
#   the toml file at the bottom of this script.
##############################################################################
##############################################################################
# NOTE: I should only save that file if each step has run without error.
##############################################################################
##############################################################################
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
            "spikenorm": None,
            "quant": None,
            "peak_calls": None,
            "epod_calls": None,
        }

    def write_ver_info(info, step, path, return_codes):
        if np.all(np.array(return_codes) == 0):
            # place version info into the preprocessing key and overwrite
            #   current file
            info[step] = VERSION 
        else:
            info[step] = "Error"
        with open(path, 'w') as f:
            toml.dump(info, f)


# define the commands that we use for each step

## The following command runs preprocessing and alignment
## It requires one argument: the config file in the working directory with detailed sample information
PR_CMD = "python {}/run_all_preprocessing.py {{}} {}".format(
    BINDIR,os.path.join(BASEDIR, conf_file)
)
AL_CMD = "python {}/run_all_alignments.py {{}} {}".format(
    BINDIR,os.path.join(BASEDIR, conf_file)
)
BS_CMD = "python {}/run_all_bootstraps.py {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)
QC_CMD = "python {}/run_all_qc.py {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)
# check whether we have paired data. If so, use one script. If not, use other.
QNORM_CMD = "python {}/quantify/qnorm_bs_files.py {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)
SPIKENORM_CMD = "python {}/quantify/spikenorm_bs_files.py {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)
QUANT_CMD = "python {}/quantify/quantify_ipod.py\
             {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)

# first just read the set of files to be acted upon
TARGETS_FILE = conf_dict_global["general"]["condition_list"]

os.chdir(BASEDIR)

instr = open(TARGETS_FILE)

all_dirs = []
all_confs = []

for line in instr:
    linearr = line.rstrip().split()
    all_dirs.append(linearr[0])
    all_confs.append(linearr[1])

# now do the first step - preprocessing and alignment

if not "preprocess" in skipsteps:
    print("Beginning preprocessing stage...")
    print("==============================================")

    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        print('')
        print('Working on sample name {}'.format(dirname))
        print('----------------------------------')
        os.chdir(os.path.join(BASEDIR, dirname))
        print(PR_CMD.format(confname))
        cp = subprocess.run(PR_CMD.format(confname), shell=True)
        retcodes.append(cp.returncode)

        print('Done with read preprocessing for sample {}'.format(dirname))
        print('____________________________________')
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "preprocessing", ver_filepath, retcodes)

if not "align" in skipsteps:
    print("Beginning alignment stage...")
    print("==============================================")

    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        print('')
        print('Working on sample name {}'.format(dirname))
        print('----------------------------------')
        os.chdir(os.path.join(BASEDIR, dirname))
        print(AL_CMD.format(confname))
        cp = subprocess.run(AL_CMD.format(confname), shell=True)
        retcodes.append(cp.returncode)

        print('Done with alignments for sample {}'.format(dirname))
        print('____________________________________')

        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "alignment", ver_filepath, retcodes)


if not "bootstrap" in skipsteps:
    print("Beginning bootstrapping stage...")
    print("==============================================")
    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        # run the quantitation commands here, making any missing directories of needed
        os.chdir(os.path.join(BASEDIR, dirname))
        print("Doing bootstrapping for sample {}".format(dirname))
        print(BS_CMD.format(confname))
        cp = subprocess.run(BS_CMD.format(confname), shell=True)
        retcodes.append(cp.returncode)
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "bootstrapping", ver_filepath, retcodes)


if not "qc" in skipsteps:
    print("Beginning QC stage...")
    print("==============================================")
    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        os.chdir(os.path.join(BASEDIR, dirname))
        print("Doing QC for sample {}".format(dirname))
        format_qc_cmd = QC_CMD.format(
            os.path.join(BASEDIR, dirname, confname)
        ) 
        print(format_qc_cmd)
        cp = subprocess.run(
            format_qc_cmd,
            shell=True
        )
        retcodes.append(cp.returncode)
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "qc", ver_filepath, retcodes)


if not "qnorm" in skipsteps:
    print("Doing quantile normalization...")
    print("==============================================")
    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        os.chdir(os.path.join(BASEDIR, dirname))
        # I need to introduce logic here to handle file names if
        #   we have paired data, since next steps assume bootstraps
        print("Doing quantile normalization for sample {}".format(dirname))
        print(QNORM_CMD.format(confname))
        cp = subprocess.run(QNORM_CMD.format(confname), shell=True)
        retcodes.append(cp.returncode)
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "qnorm", ver_filepath, retcodes)


if not "spikenorm" in skipsteps:
    print("Doing spike-in normalization...")
    print("==============================================")
    retcodes = []
    for dirname,confname in zip(all_dirs, all_confs):
        os.chdir(os.path.join(BASEDIR, dirname))
        print("Doing spike-in normalizations for sample {}".format(dirname))
        print(SPIKENORM_CMD.format(confname))
        cp = subprocess.run(SPIKENORM_CMD.format(confname))
        retcodes.append(cp.returncode)
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "spikenorm", ver_filepath, retcodes)



if not "quant" in skipsteps:
    print("\nEntering quantitation step")
    print("==============================================")
    retcodes = []

    for dirname,confname in zip(all_dirs, all_confs):
        print('')
        print('Working on sample name {}'.format(dirname))
        print('----------------------------------')

        # run the quantitation commands here, making any missing directories of needed
        os.chdir(os.path.join(BASEDIR, dirname))
        print("== running quant cmd")
        print("== running quant cmd", file=sys.stderr)
        print(QUANT_CMD.format(confname))
        cp = subprocess.run(QUANT_CMD.format(confname), shell=True)
        retcodes.append(cp.returncode)

        print('Done processing sample {}'.format(dirname))
        print('____________________________________')

        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "quant", ver_filepath, retcodes)

print("==============================================")
print("FINISHED WITH ALL SAMPLES")
