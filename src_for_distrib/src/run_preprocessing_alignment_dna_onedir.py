#!/usr/bin/python

# script for doing my standardized preprocessing and alignment on DNA files
# we do standardized quality trimming, adapter removal, and then alignment to the specified genome
# we expect a fairly standardized directory layout:
#  There should be  'aligned' subdirectory to contain the aligned data
#  If there is not already an existing 'raw' directory, we link the general->rawdir file here
#  In addition, the top level directory should have a space-delimited file called 'seq_manifest.txt' 
#  This file needs to have a series of input prefix/output prefix pairs, and will be used to identify which runs to look at
#  There are two command line arguments: 
#    the name of the sample manifest file to be read to find sample information
#    the path to the top-level config file containing genome information

import sys
import subprocess
import tempfile
import os
import toml

# parse the top-level config and get needed general information
conf_dict = toml.load(sys.argv[2])

## Set up defined constants that should be universal
proc_opts = conf_dict["processing"]
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

aln_opts = conf_dict["alignment"]
ALDIR = aln_opts["aligned_direc"]
NPROC = aln_opts["align_threads"]
MIN_FRAG_LEN = aln_opts["min_fragment_length"]
MAX_FRAG_LEN = aln_opts["max_fragment_length"]

SEQ_DB = conf_dict["genome"]["genome_base"]
# base directory for all code for ipod analysis
BINDIR = conf_dict["general"]["bindir"]

# set up the needed directories if they are not already present
if not(os.path.isdir(PROCDIR)):
    os.mkdir(PROCDIR)
if not(os.path.isdir(ALDIR)):
    os.mkdir(ALDIR)

if not(os.path.islink('raw')):
    os.symlink(conf_dict["general"]["rawdir"], "raw")

# the actual input to this program should be a single space-delmited file
# each line should have, in order:
# Freads Rreads adapseq phredbase outprefix

#now, we go parse that file
all_samples = []

instr = open(sys.argv[1]) # this is the read_manifest.txt file

for line in instr:

    if line[0] == '#':
        continue

    linearr = line.rstrip().split()
    this_samp = dict()
    this_samp["ffile"] = linearr[0]
    this_samp["rfile"] = linearr[1]
    this_samp["adapseq"] = linearr[2]
    this_samp["phredbase"] = int(linearr[3])
    this_samp["outprefix"] = linearr[4]
    all_samples.append(this_samp)

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
    trim_cmd = "TrimmomaticPE -threads {} -phred{} \
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

# now we can actually run the samples

for samp in all_samples:
    print("running on {}".format(samp))
    # NOTE: could rewrite for multiprocessing.
    preprocess_gz_file(samp)

for samp in all_samples:
    run_bowtie(samp["outprefix"],samp["phredbase"],db=SEQ_DB)
    postprocess_bowtie(samp["outprefix"])
