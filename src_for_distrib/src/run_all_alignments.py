#!/usr/bin/python

# run all preprocessing and alignments for all samples specified in a given config file
# this is intended to be run from the top level directory for those samples
# the config file should indicate in the [general]->sample_types option what types of samples (inp, chip, etc) are present; there
#  should already be a directory for each such sample type.
# in addition, each of those directories should have a properly filled out read_manifest.txt, and must have raw/, aligned/, and bootstrap/ directories 
# already set up
# the command line arguments are:
#  -a config file within the directory of interest indicating where all of the samples are
#  -the top level config file that specifies general information like genome size, genome location, etc.

import sys
import subprocess
import toml
import os

# parse the top level config file and get needed information

BASE_DIR = os.getcwd()

conf_file = sys.argv[1]
conf_file_global = sys.argv[2]
conf_dict = toml.load(conf_file)
conf_dict_global = toml.load(conf_file_global)

ALIGN_CMD = "python {}/run_preprocessing_alignment_dna_onedir.py\
             read_manifest.txt {}".format(
    conf_dict_global["general"]["bindir"], conf_file_global
)

n_errors = 0
for sample_type in conf_dict["general"]["sample_types"]:
    print("Working on {} samples...".format(sample_type))
    try:
        os.chdir(conf_dict[sample_type]["directory"])
        # Now we're in the directory containing the fastq files
        #   and read_manifest.txt file.
        subprocess.run(ALIGN_CMD, shell=True, check=True)
        # after ALIGN_CMD completes, go back to BASE_DIR
        os.chdir(BASE_DIR)
    except:
        n_errors += 1
        print(
            "Warning: Encountered an error processing {} samples".format(
                sample_type
            )
        )


print(
    "Finished with preprocessing and alignment. Encountered {} errors".format(
        n_errors
    )
)
