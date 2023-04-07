#!/usr/bin/python

# actually run all of the needed commands in this directory

import sys
import os 
import shutil
import toml
import argparse
import glob

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
    '--cleansteps',
    help="comma-separated list of steps, the files for which to remove\
         (can be any of preprocess,align,bootstrap)",
    default=None
)

args = parser.parse_args()
conf_file = args.main_conf
# cleansteps defaults to empty set
if args.cleansteps is None:
    cleansteps = set()
else:
    cleansteps = set(args.cleansteps.split(","))

steps = ['preprocess','align','bootstrap']
for step in cleansteps:
    if step not in steps:
        raise NoSuchStepException(step, steps)

# parse the top level config file to get some needed information
conf_dict_global = toml.load(conf_file)

BASEDIR = conf_dict_global["general"]["basedir"]
BINDIR = conf_dict_global["general"]["bindir"]

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
if "preprocess" in cleansteps:
    print("Cleaning up files created during the preprocessing stage...")
    print("===========================================================")

    for dirname,confname in zip(all_dirs, all_confs):
        os.chdir(os.path.join(BASEDIR, dirname))
        cond_conf = toml.load(confname)
        samptypes = cond_conf["general"]["sample_types"]
        for samptype in samptypes:
            os.chdir(os.path.join(
                samptype,conf_dict_global["processing"]["processed_direc"]
            ))
            for fqfile in glob.glob("*.fq.gz"):
                os.remove(fqfile)

        os.chdir(BASEDIR)

if "align" in cleansteps:
    print("Cleaning up files created during the alignment stage...")
    print("===========================================================")

    for dirname,confname in zip(all_dirs, all_confs):
        os.chdir(os.path.join(BASEDIR, dirname))
        cond_conf = toml.load(confname)
        samptypes = cond_conf["general"]["sample_types"]
        for samptype in samptypes:
            os.chdir(os.path.join(
                samptype,conf_dict_global["alignment"]["aligned_direc"]
            ))
            for fname in glob.glob("*.bam*"):
                os.remove(fname)

        os.chdir(BASEDIR)

