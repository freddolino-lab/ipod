#!/usr/bin/python

# actually run all of the needed commands in this directory

import subprocess
import sys
import os 
import toml
import argparse
import numpy as np

class NoSuchStepException(Exception):
    def __init__(self, step, steps):
        self.message = f"ERROR: {step} is not an allowable "\
            f"step to skip. Allowed steps are "\
            f"{' '.join(steps)}."
        super().__init__(self.message)

# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    'main_conf',
    help="main configuration file defining work to be done"
)
parser.add_argument(
    '--skipsteps',
    help="comma-separated list of steps to skip\
         (can be any of preprocess,align,bootstrap,qc,qnorm,spikenorm,quant)",
    default=None
)

args = parser.parse_args()
conf_file = args.main_conf
# skipsteps defaults to empty set
if args.skipsteps is None:
    skipsteps = set()
else:
    skipsteps = set(args.skipsteps.split(","))

steps = ['umi','preprocess','align','bootstrap','qc','qnorm','spikenorm','quant']
for step in skipsteps:
    if step not in steps:
        raise NoSuchStepException(step, steps)

# parse the top level config file to get some needed information
conf_dict_global = toml.load(conf_file)

BASEDIR = conf_dict_global["general"]["basedir"]
BINDIR = conf_dict_global["general"]["bindir"]
SPIKE = False
SPIKE_CHR = conf_dict_global["genome"]["spike_in_name"]
if SPIKE_CHR != "None":
    SPIKE = True

# if we've run this driver from within our singularity container,
#   then the IPOD_VER environment veriable with exist.
# In that case, read in the toml file (if it exists) denoting which versions of
#   the container were used for which steps in the past. If the file doesn't 
#   exist, create a dictionary populated with correct information, and save
#   the toml file at the bottom of this script.
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
            sys.exit(f"Error in {step}. Exiting now. Check logs.")
        with open(path, 'w') as f:
            toml.dump(info, f)


# define the commands that we use for each step
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

        if cp.returncode != 0:
            sys.exit(f"Error in preprocessing for {dirname}. Check logs. Exiting now.")

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

        if cp.returncode != 0:
            sys.exit(f"Error in alignment for {dirname}. Check logs. Exiting now.")

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

        if cp.returncode != 0:
            sys.exit(f"Error in bootstrapping for {dirname}. Check logs. Exiting now.")

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

        if cp.returncode != 0:
            sys.exit(f"Error in qc for {dirname}. Check logs. Exiting now.")

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
        print("Doing quantile normalization for sample {}".format(dirname))
        print(QNORM_CMD.format(confname))
        cp = subprocess.run(QNORM_CMD.format(confname), shell=True)

        if cp.returncode != 0:
            sys.exit(f"Error in quantile normalization for {dirname}. Check logs. Exiting now.")

        retcodes.append(cp.returncode)
        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "qnorm", ver_filepath, retcodes)


if not "spikenorm" in skipsteps:
    if SPIKE:
        print("Doing spike-in normalization...")
        print("==============================================")
        retcodes = []
        for dirname,confname in zip(all_dirs, all_confs):
            os.chdir(os.path.join(BASEDIR, dirname))
            print("Doing spike-in normalizations for sample {}".format(dirname))
            print(SPIKENORM_CMD.format(confname))
            cp = subprocess.run(SPIKENORM_CMD.format(confname), shell=True)

            if cp.returncode != 0:
                sys.exit(f"Error in spike-in normalization for {dirname}. Check logs. Exiting now.")

            retcodes.append(cp.returncode)
            os.chdir(BASEDIR)

        if "IPOD_VER" in os.environ:
            write_ver_info(ver_info, "spikenorm", ver_filepath, retcodes)
    else:
        print("The name of the spike-in normalization contig was set to 'None' in your main configuration file. Skipping spike-in normalization.")


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

        if cp.returncode != 0:
            sys.exit(f"Error in quan step for {dirname}. Check logs. Exiting now.")

        retcodes.append(cp.returncode)

        print('Done processing sample {}'.format(dirname))
        print('____________________________________')

        os.chdir(BASEDIR)

    if "IPOD_VER" in os.environ:
        write_ver_info(ver_info, "quant", ver_filepath, retcodes)

print("==============================================")
print("FINISHED WITH ALL SAMPLES")
