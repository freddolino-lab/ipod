#!/usr/bin/python

# Run all bootstrapping needed for the samples in a given config file
# This is intended to be run from the top level directory for those samples
# The config file should indicate in the [general]->sample_types option what
#   types of samples (inp, chip, etc) are present; there should already be
#   a directory for each such sample type.
# In addition, each of those directories should have a properly filled out
#   read_manifest.txt, and must have raw/, aligned/, and bootstrap/ directories 
#   already set up
# We then do a standard set of bootstrapping where we generate both occupancy
#   traces for the original data, and a specified number
#   of bootstrap replicates.
# Command line arguments are the sample-specific and experiment-wide config
#   files, in order.

import sys
import subprocess
import toml
import os
import multiprocessing
import tempfile
import pickle
import h5py
import numpy as np
import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, 'utils')
sys.path.insert(0, utils_path)

import hdf_utils

conf_file = sys.argv[1]
conf_dict = toml.load(conf_file)
conf_file_global = sys.argv[2]
conf_dict_global = toml.load(conf_file_global)

BINDIR = conf_dict_global["general"]["bindir"]
BSDIR = conf_dict_global["bootstrap"]["bootstrap_direc"]
BS_SAM_THREADS = conf_dict_global["bootstrap"]["samtools_threads"]
BS_BOOT_THREADS = conf_dict_global["bootstrap"]["bootstrap_threads"]

# GENOME_SIZE=conf_dict_global["genome"]["genome_length"]
BS_RESOLUTION = conf_dict_global["genome"]["resolution"]
BS_NUMSAMP = conf_dict_global["bootstrap"]["bootstrap_samples"]
RETAIN_FLAGS = conf_dict_global["bootstrap"]["aln_retain_flags"]
REJECT_FLAGS = conf_dict_global["bootstrap"]["aln_reject_flags"]
MAPQ = conf_dict_global["bootstrap"]["aln_mapq_filter"]

BS_SUFFIX = conf_dict_global["bootstrap"]["bs_suffix"]

# get genome size from bowtie2 index
# NOTE: somthing like this will have to be done do get multi-contig info
inspect_cmd = 'bowtie2-inspect -s {} | grep ^Sequence | cut -f 2,3'.format(
    conf_dict_global["genome"]["genome_base"]
)

CONTIG_SIZES = subprocess.check_output(inspect_cmd, shell=True)
CONTIG_SIZES = CONTIG_SIZES.decode().strip().split('\n')
# Sorting inplace here will ensure indices for lookup table are
#   sorted according to alpha sorting of contig names.
CONTIG_SIZES.sort()

# PARSE_CMD will get hdf_file and input samfile substituted in later.
PARSE_CMD = "python {}/bootstrapping/bootstrap_sam_file.py\
                 --hdf_file {{}} --global_conf_file {} parse\
                 --paired {{}}".format(BINDIR, conf_file_global)
# SAMPLE_CMD and ORIG_CMD will get hdf_file substituted later
SAMPLE_CMD = "python {}/bootstrapping/bootstrap_sam_file.py\
               --hdf_file {{}} --global_conf_file {} sample --num_samples {}\
               --resolution {}".format(
    BINDIR, conf_file_global, BS_NUMSAMP, BS_RESOLUTION
)
ORIG_CMD = "python {}/bootstrapping/bootstrap_sam_file.py\
             --hdf_file {{}} --global_conf_file {} sample --identity\
             --resolution {}".format(
    BINDIR, conf_file_global, BS_RESOLUTION
)
n_errors = 0

# here we define a helper function to do the bootstrap part ONLY
def preprocess_bootstrap( bamfile, hdf_name, nthreads ):
    '''Do the preprocessing necessary to set up a bam file for 
    bootstrapping. This does the initial setup and parse stages, 
    but not the actual bootstrap resampling.

    Args:
    ----
    bamfile : str
        The path to the sorted bamfile
    outprefix : str
        Prefix (sample name) to prepend to outputs
    nthreads : int
        Number of threads to run samtools with simultaneously

    Returns:
    -------
    None
        Generates a npy file containing an Nx2 array, where N
        is the number of reads passing the filters defined
        in the "bootstrap" section of the main configuration file.
        For each row (read) the [start,end] positions are recorded
        in the npy file.
    '''

    # first get a sorted sam file as needed by the parse command
    tmp_file = tempfile.NamedTemporaryFile()
    print((PARSE_CMD.format(hdf_name, tmp_file.name)))
    filter_cmd = "samtools sort -n -@ {} {} \
                    | samtools view -f {} -F {} -q {} - \
                    | sed s/\"#0\/4\t\"/\"#0\t\"/g > {}".format(
        nthreads, bamfile, RETAIN_FLAGS, REJECT_FLAGS, MAPQ, tmp_file.name
    )
    print("\n{}\n".format(filter_cmd))
    subprocess.call(filter_cmd, shell=True)

    # now run the parse step
    subprocess.call(
        PARSE_CMD.format(hdf_name, tmp_file.name),
        shell=True
    )
    # and clean up
    tmp_file.close()


def do_bootstrap( hdf_name ):
    '''Assuming that parsing has already completed for a specified bam file,
       now we actually run the needed resampling.

    Args:
    -----
    outprefix : str
        Prefix (sample name) to prepend to outputs
        
    Returns:
    --------
    None
        Generates two npy files.
        The first npy file contains a numpy array of shape (G,B), where G is 
        the genome length and B is the number of bootstrap replicates.
        Values of the first array are the bootstrap-sampled coverages
        for each position of the reference genome (rows), for each 
        replicate (columns).
        The second npy file contains the observed coverages for each
        position of the reference genome.
    '''

    # do the actual bootstrapping 
    print(SAMPLE_CMD.format(hdf_name))
    subprocess.call(
        SAMPLE_CMD.format(hdf_name),
        shell=True
    )

    # also generate a file with the original occupancies
    print(ORIG_CMD.format(hdf_name))
    subprocess.call(
        ORIG_CMD.format(hdf_name),
        shell=True
    )

# final driver function to actually do all of the work needed for each
#   individual sample type
def bs_sample_type(samptype, sampledir, ctg_lut, n_errors):
    '''Main driver function to run all bootstrapping for a
    specified sample type.

    Args:
    -----
    samptype : str
        Can be any of {inp,chip,nitr,ipod}. Denotes whether this sample is from
        input DNA, ipod extraction, chip-seq, or nirocellulose-binding ipod.
    sampledir : str
        Directory containing the samples to be bootstrapped.
    ctg_lut : dict
        Dictionary with information on contig ids and sizes to be used in
        setting up hdf5 files to store multi-contig parsers, coverages, and
        bootstraps.
    n_errors : int
        Iterator to count the number of errors the script encountered during
        this run.

    Returns:
    --------
    None
        This function wraps the parsing and resampling steps into one function.
        It also appends the resulting processes to the list of threads.
    '''
    print("Working on {} samples...".format(samptype))

    sample_prefixes = conf_dict[samptype]["sample_names"]
    bam_files = [
        os.path.join(
            sampledir,
            conf_dict_global["alignment"]["aligned_direc"],
            pref + "_bowtie2_sorted.bam" )
        for pref in sample_prefixes
    ]
    out_names = [
        os.path.join( sampledir, BSDIR, pref )
        for pref in sample_prefixes
    ]

    # make the bootstrap directory if it does not already exist
    if not(os.path.isdir(os.path.join(sampledir, BSDIR))):
        os.mkdir(os.path.join(sampledir, BSDIR))

    for bamname, outpref in zip(bam_files, out_names):

        print("Preprocessing {}".format(bamname))

        hdf_name = outpref + ".hdf5"
        hdf_utils.set_up_hdf_file(hdf_name, ctg_lut, BS_RESOLUTION)

        preprocess_bootstrap(bamname, hdf_name, BS_SAM_THREADS)

        print("Background processing bootstrap for {}".format(bamname))
        all_thr.append(bs_pool.apply_async( do_bootstrap, [hdf_name] ))

# Now set up the pool that we will use for most of the work.
bs_pool = multiprocessing.Pool(BS_BOOT_THREADS)
all_thr = [] # this is a list that will contain all of the result objects

# set up lookup table for contig indices to their id and length
ctg_lut = {}
for ctg_idx,ctg in enumerate(CONTIG_SIZES):

    ctg_id = ctg.split(' ')[0]
    ctg_len = int(ctg.split('\t')[-1])
    ctg_lut[ctg_idx] = {"id":ctg_id, "length":ctg_len}

print("Beginning work on bootstrapping")
n_errors = 0
for sample_type in conf_dict["general"]["sample_types"]:
    bs_sample_type(
        sample_type,
        conf_dict[sample_type]["directory"],
        ctg_lut,
        n_errors,
    )

# now collect all of the threads
bs_pool.close()
bs_pool.join()

n_pool_err = 0
for res in all_thr:
    if not res.successful():
        n_pool_err += 1

print("Finished with bootstrapping. \
    Encountered {} errors with preprocessing and \
    {} errors during the bootstrap portion".format(
        n_errors, n_pool_err
))
