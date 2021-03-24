#!/usr/bin/python

# actually run all of the needed commands in this directory

import subprocess
import sys
import os 
import toml
import argparse

# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    'main_conf',
    help="main configuration file defining work to be done"
)
parser.add_argument(
    '--skipsteps',
    help="comma-separated list of steps to skip\
         (can be any of align,bootstrap,qc,qnorm,quant)",
    default=None
)

args = parser.parse_args()
conf_file = args.main_conf

if args.skipsteps is None:
    skipsteps = set()
else:
    skipsteps = set(args.skipsteps.split(","))

# parse the top level config file to get some needed information
conf_main = toml.load(conf_file)

BASEDIR = conf_main["general"]["basedir"]
BINDIR = conf_main["general"]["bindir"]

# define the commands that we use for each step

## The following command runs preprocessing and alignment
## It requires one argument: the config file in the working directory with detailed sample information
PR_AL_CMD = "python {}/run_all_alignments.py {{}} {}".format(
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
QUANT_CMD = "python {}/quantify/quantify_ipod.py\
             {{}} {}".format(
    BINDIR, os.path.join(BASEDIR, conf_file)
)

# first just read the set of files to be acted upon
TARGETS_FILE = conf_main["general"]["condition_list"]

instr = open(TARGETS_FILE)

all_dirs = []
all_confs = []

for line in instr:
    linearr = line.rstrip().split()
    all_dirs.append(linearr[0])
    all_confs.append(linearr[1])

# now do the first step - preprocessing and alignment

if not "align" in skipsteps:
    print("Beginning preprocessing and alignment stage...")
    print("==============================================")

    for dirname,confname in zip(all_dirs, all_confs):
        print('')
        print('Working on sample name {}'.format(dirname))
        print('----------------------------------')

        os.chdir(os.path.join(BASEDIR, dirname))
        print(PR_AL_CMD.format(confname))
        subprocess.run(PR_AL_CMD.format(confname), shell=True)

        print('Done with alignments for sample {}'.format(dirname))
        print('____________________________________')

        os.chdir(BASEDIR)

print("Beginning bootstrapping and QC stage...")
print("==============================================")


for dirname,confname in zip(all_dirs, all_confs):
    print('')
    print('Working on sample name {}'.format(dirname))
    print('----------------------------------')

    # run the quantitation commands here, making any missing directories of needed
    os.chdir(os.path.join(BASEDIR, dirname))
    if not "bootstrap" in skipsteps:
        print("Doing bootstrapping for sample {}".format(dirname))
        print(BS_CMD.format(confname))
        subprocess.run(BS_CMD.format(confname), shell=True)
    if not "qc" in skipsteps:
        print("Doing qc")
        print(QC_CMD.format(confname))
        subprocess.run(
            QC_CMD.format(confname),
            shell=True
        )
    if not "qnorm" in skipsteps:
        # I need to introduce logic here to handle file names if
        #   we have paired data, since next steps assume bootstraps
        print("Doing quantile normalization for sample {}".format(dirname))
        print(QNORM_CMD.format(confname))
        subprocess.run(QNORM_CMD.format(confname), shell=True)

    print('____________________________________')

    os.chdir(BASEDIR)


if not "quant" in skipsteps:
    print("\nEntering quantitation step")
    print("==============================================")

    for dirname,confname in zip(all_dirs, all_confs):
        print('')
        print('Working on sample name {}'.format(dirname))
        print('----------------------------------')

        # run the quantitation commands here, making any missing directories of needed
        os.chdir(os.path.join(BASEDIR, dirname))
        print("== running quant cmd")
        print("== running quant cmd", file=sys.stderr)
        print(QUANT_CMD.format(confname))
        subprocess.run(QUANT_CMD.format(confname), shell=True)

        print('Done processing sample {}'.format(dirname))
        print('____________________________________')

        os.chdir(BASEDIR)

print("==============================================")
print("FINISHED WITH ALL SAMPLES")
