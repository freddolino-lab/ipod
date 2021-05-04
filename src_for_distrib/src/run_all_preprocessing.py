#!/usr/bin/python

# script for doing my standardized preprocessing on DNA files
# we do standardized quality trimming and adapter removal

import sys
import subprocess
import tempfile
import os
import toml

# parse the top-level config and get needed general information
conf_dict_global = toml.load(sys.argv[2])

## Set up defined constants that should be universal
proc_opts = conf_dict_global["processing"]
PROCDIR = proc_opts["processed_direc"]
MAX_ADAPT_N = proc_opts["adapt_max_n"]
TRAILING_JUNK_LEN = proc_opts["trim_trailing_junk_length"]
SLIDE_WIN_LEN = proc_opts["trim_sliding_window_length"]
SLIDE_WIN_QUAL = proc_opts["trim_sliding_window_qual"]
PROCESSED_READ_MINLEN = proc_opts["min_processed_readlen"]
F_READ_SUFFIX = proc_opts["f_paired_read_file_suffix"]
R_READ_SUFFIX = proc_opts["r_paired_read_file_suffix"]
F_UP_READ_SUFFIX = proc_opts["f_unpaired_read_file_suffix"]
R_UP_READ_SUFFIX = proc_opts["r_unpaired_read_file_suffix"]

aln_opts = conf_dict_global["alignment"]
NPROC = aln_opts["align_threads"]
MIN_FRAG_LEN = aln_opts["min_fragment_length"]
MAX_FRAG_LEN = aln_opts["max_fragment_length"]

# base directory for all code for ipod analysis
BINDIR = conf_dict_global["general"]["bindir"]
RAWDIR = conf_dict_global["general"]["rawdir"]
STARTDIR = os.getcwd()

# if the rawdir option is not "None" (note a string, not a None object)
#   then determine whether it's already a symlink. If it's not a symlink
#   already, create the symlink within the data directory
if RAWDIR != "None":
    if not(os.path.islink('raw')):
        os.symlink(conf_dict_global["general"]["rawdir"], "raw")

# define some functions that will be used in the rest of the script
def preprocess_gz_file(samp):
    '''Do some initial preprocessing of a gz file,
    including trimming and quality score filtering.
    Don't use this function outside of this script;
    it depends on information defined above.

    Args:
    -----
    samp : dict
        This dict object has all the data listed
        higher up in this file

    Returns:
    --------
    None
        Creates fq.gz files containing processed reads
        Two files for the paired reads (fr, rev)
        Two files for the unpaired reads (fr, rev)
    '''

    infile_1 = samp["ffile"]
    infile_2 = samp["rfile"]
    outprefix = samp["outprefix"]
    PHRED_BASE = samp["phredbase"]
    ADAP_SEQ = samp["adapseq"]

    if infile_1[-3:] == ".gz":
        DCPROG = 'zcat'
    elif infile_1[-4:] == ".bz2":
        DCPROG = 'bzcat'
    else:
          raise("Could not determine the decompression program to use")

    if DCPROG == "bzcat":
        in1 = tempfile.NamedTemporaryFile(suffix='.fastq')
        in2 = tempfile.NamedTemporaryFile(suffix='.fastq')
        infile_fwd = in1.name
        infile_rev = in2.name
        cmd1="{} {} > {}".format(DCPROG, infile_1, infile_fwd)
        cmd2="{} {} > {}".format(DCPROG, infile_2, infile_rev)
        subprocess.call(cmd1,shell=True)
        subprocess.call(cmd2,shell=True)

    else:
        infile_fwd = infile_1
        infile_rev = infile_2

    cutfile_fwd = os.path.join(PROCDIR, outprefix+"_fwd_cutadap.fq.gz")
    cutfile_rev = os.path.join(PROCDIR, outprefix+"_rev_cutadap.fq.gz")

    #  do some quality trimming and write a processed file
    # first clip the adapter sequences
    cutadapt_cmd = "cutadapt --quality-base={} \
                             -a {} -A {} -n {} --match-read-wildcards \
                             -o {} -p {} {} {} \
                             > {}_cutadapt.log 2> {}_cutadapt.err".format(
        PHRED_BASE, ADAP_SEQ, ADAP_SEQ, MAX_ADAPT_N,
        cutfile_fwd, cutfile_rev, infile_fwd, infile_rev,
        outprefix, outprefix
    )
    print("\n{}\n".format(cutadapt_cmd))
    subprocess.call(cutadapt_cmd,shell=True)

    # Next do quality trimming -- trim true crap from the 3' end,
    #   and then look for a  sliding window of 4 bp with qualities above 15.
    # We avoid doing 5' end trimming so we don't move where the true end
    #  of the read is.
    # Drop the read if we have less than 10 bases after this is done
    trim_fwd_paired = os.path.join( PROCDIR, outprefix + F_READ_SUFFIX )
    trim_fwd_unpaired = os.path.join( PROCDIR, outprefix + F_UP_READ_SUFFIX )
    trim_rev_paired = os.path.join( PROCDIR, outprefix + R_READ_SUFFIX )
    trim_rev_unpaired = os.path.join( PROCDIR, outprefix + R_UP_READ_SUFFIX )
    trim_cmd = "trimmomatic PE -threads {} -phred{} \
                    {} {} {} {} {} {} \
                    TRAILING:{} SLIDINGWINDOW:{}:{} \
                    MINLEN:{} 2> {}_trimmomatic.err".format(
        NPROC, PHRED_BASE, cutfile_fwd, cutfile_rev,
        trim_fwd_paired, trim_fwd_unpaired, trim_rev_paired,
        trim_rev_unpaired, TRAILING_JUNK_LEN, SLIDE_WIN_LEN,
        SLIDE_WIN_QUAL, PROCESSED_READ_MINLEN, outprefix
    )
    print("\n{}\n".format(trim_cmd))
    subprocess.call(trim_cmd,shell=True)

    if DCPROG == "bzcat":
        in1.close()
        in2.close()


# the actual input to this program should be a single space-delmited file
# each line should have, in order:
# Freads Rreads adapseq phredbase outprefix

conf_dict = toml.load(sys.argv[1]) # this is the condition-level conf file
samp_types = conf_dict["general"]["sample_types"]

for samp_type in samp_types:

    print("Running preprocessing for {} sample.".format(samp_type))

    # gather this sample type's info from conf file
    samp_info = conf_dict[samp_type]
    freads = samp_info["R1_raw_files"]
    rreads = samp_info["R2_raw_files"]
    adapts = samp_info["adapter_seqs"]
    rep_names = samp_info["sample_names"]
    sample_direc = samp_info["directory"]

    # move into this sample's directory
    os.chdir(sample_direc)

    # set up the needed directories if they are not already present
    if not(os.path.isdir(PROCDIR)):
        os.mkdir(PROCDIR)

    # loop over each replicate's information and do preprocessing
    for i in range(len(rep_names)):
        samp_dict = {
            "ffile": freads[i],
            "rfile": rreads[i],
            "adapseq": adapts[i],
            "phredbase": conf_dict_global["general"]["phredbase"],
            "outprefix": rep_names[i],
        }
        preprocess_gz_file(samp_dict)
    
    # when finished with this sample's preprocessing, move back up
    #  to condition-level directory
    os.chdir(STARTDIR)

