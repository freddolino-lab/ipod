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

def clear_dir(conf_dir, confname, clearstep, cleardir, fnamesearch):
    conf_file = os.path.join(conf_dir, confname)
    cond_conf = toml.load(conf_file)
    samptypes = cond_conf["general"]["sample_types"]
    for samptype in samptypes:
        samp_dir = cond_conf[samptype]["directory"]
        this_dir = os.path.join(
            conf_dir, samp_dir, conf_dict_global[clearstep][cleardir]
        )
        search_term = os.path.join(this_dir, fnamesearch)
        if os.path.isdir(this_dir):
            print(f"\nLooking in {this_dir} for files ending in {fnamesearch}.")
            to_remove = glob.glob(search_term)
            #print(f"Matched files marked for removal:\n{to_remove}")
            if not to_remove:
                print(f"No files found matching {search_term}. Moving on.")
            for fname in to_remove:
                print(f"Removing file {fname}")
                os.remove(fname)
        else:
            print(f"Currently in {os.getcwd()}.")
            print(f"No directory {this_dir}. Skipping it.")
            continue



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
        this_dir = os.path.join(BASEDIR, dirname)

        #os.chdir(this_dir)
        clear_dir(this_dir, confname, "processing", "processed_direc", "*.fq.gz")

        #os.chdir(BASEDIR)

if "align" in cleansteps:
    print("Cleaning up files created during the alignment stage...")
    print("===========================================================")

    for dirname,confname in zip(all_dirs, all_confs):
        this_dir = os.path.join(BASEDIR, dirname)
        #os.chdir(os.path.join(BASEDIR, dirname))
        clear_dir(this_dir, confname, "alignment", "aligned_direc", "*.bam*")
        #os.chdir(BASEDIR)

if "bootstrap" in cleansteps:
    print("Cleaning up files created during the read bootstrapping stage...")
    print("===========================================================")

    for dirname,confname in zip(all_dirs, all_confs):
        this_dir = os.path.join(BASEDIR, dirname)
        #os.chdir(os.path.join(BASEDIR, dirname))
        clear_dir(this_dir, confname, "bootstrap", "bootstrap_direc", "*.hdf5")
        #os.chdir(BASEDIR)



