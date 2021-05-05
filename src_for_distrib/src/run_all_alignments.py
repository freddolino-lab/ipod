#!/usr/bin/python

# script for doing my standardized preprocessing and alignment on DNA files
# we do standardized quality trimming, adapter removal, and then alignment to the specified genome

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
ALDIR = aln_opts["aligned_direc"]
NPROC = aln_opts["align_threads"]
MIN_FRAG_LEN = aln_opts["min_fragment_length"]
MAX_FRAG_LEN = aln_opts["max_fragment_length"]

SEQ_DB = conf_dict_global["genome"]["genome_base"]
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
def run_bowtie(prefix, phredbase, db=SEQ_DB):
    '''Run alignment using bowtie2

    Args:
    -----
    prefix : str
        Characters to prepend to files generated during alignment
    phredbase : int
        Quality score encoding
    db : str
        Path to bowtie2 database to serve as reference for alignment

    Returns:
    --------
    None
        Runs alignment and generates sam file
    '''
    fwd = os.path.join(PROCDIR,prefix+F_READ_SUFFIX)
    rev = os.path.join(PROCDIR,prefix+R_READ_SUFFIX)
    fwd_unpaired = os.path.join( PROCDIR, prefix + F_UP_READ_SUFFIX )
    rev_unpaired = os.path.join( PROCDIR, prefix + R_UP_READ_SUFFIX )
    samout = os.path.join(ALDIR,prefix+"_bowtie2.sam")
    cmdline = 'bowtie2 -x {} -1 {} -2 {} \
                       -U {},{} -S {} \
                       -q --end-to-end --very-sensitive \
                       -p {} --no-unal --phred{} \
                       --fr -I {} -X {} \
                       --dovetail > {}_bowtie2.log 2> {}_bowtie2.err'.format(
        db, fwd, rev, fwd_unpaired,
        rev_unpaired, samout, NPROC, phredbase,
        MIN_FRAG_LEN, MAX_FRAG_LEN, prefix, prefix
    )

    print("\n{}\n".format(cmdline))
    subprocess.call(cmdline, shell=True)
  

def postprocess_bowtie(prefix):
    '''Convert sam file to bam and sort the resulting bam file

    Args:
    -----
    prefix : str
        Prefix to prepend to generated files

    Returns:
    --------
    None
        When finished we're left with a sorted, unfiltered bam file and
        its associated index (bai) file.
    '''
    samname = os.path.join(ALDIR, prefix+"_bowtie2.sam")
    bamname_un = os.path.join(ALDIR, prefix+"_unsorted.bam")
    bamname = os.path.join(ALDIR, prefix+"_bowtie2_sorted.bam")
    cmdline1 = "samtools view -bS {} > {}".format(samname, bamname_un)
    cmdline2 = "samtools sort -o {} {}".format(bamname, bamname_un)
    cmdline3 = "samtools index {}".format(bamname)
    print("\n{}\n".format(cmdline1))
    print("\n{}\n".format(cmdline2))
    print("\n{}\n".format(cmdline3))
    retcode1 = subprocess.call(cmdline1,shell=True)
    retcode2 = subprocess.call(cmdline2,shell=True)
    retcode3 = subprocess.call(cmdline3,shell=True)
    if retcode1 == 0 and retcode2 == 0 and retcode3 == 0:

        print("Safe to remove samfile {}".format(prefix))
        os.remove(samname)
        os.remove(bamname_un)
    else:
        print("*** Encountered an error while postprocessing {}".format(prefix))

conf_dict = toml.load(sys.argv[1]) # this is the condition-level conf file
samp_types = conf_dict["general"]["sample_types"]

for samp_type in samp_types:

    print("Running alignments for {} sample.".format(samp_type))
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
    if not(os.path.isdir(ALDIR)):
        os.mkdir(ALDIR)
   
    # loop over each replicate's information and do preprocessing
    for i in range(len(rep_names)):
        samp_dict = {
            "ffile": freads[i],
            "rfile": rreads[i],
            "adapseq": adapts[i],
            "phredbase": conf_dict_global["general"]["phredbase"],
            "outprefix": rep_names[i],
        }
        run_bowtie(
            samp_dict["outprefix"],
            samp_dict["phredbase"],
            db=SEQ_DB,
        )
        postprocess_bowtie(samp_dict["outprefix"])

    # when finished with this sample's alignments, move back up
    #  to condition-level directory
    os.chdir(STARTDIR)

