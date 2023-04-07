#!/usr/bin/python

# script for doing my standardized preprocessing and alignment on DNA files
# we do standardized quality trimming, adapter removal, and then alignment to the specified genome

import sys
import subprocess
import tempfile
import os
import toml

class NoRefException(Exception):
    def __init__(self, db):
        self.message = f"ERROR: bowtie2 was unable to locate your reference at {db}. "\
            f"Did you mount the location containing your reference, "\
            f"and does your main configuration file point to it correctly "\
            f"within your container?"
        super().__init__(self.message)

# parse the top-level config and get needed general information
conf_dict_global = toml.load(sys.argv[2])

## Set up defined constants that should be universal
proc_opts = conf_dict_global["processing"]
PE = conf_dict_global["general"]["paired_reads"]
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
# set to false by default, and switch only if handle_umi is there, and is set to true.
# this provides backward-compatability
UMI = False
if "handle_umi" in proc_opts:
    UMI = proc_opts["handle_umi"]
    UMI_METHOD = conf_dict_global["umi"]["method"]

aln_opts = conf_dict_global["alignment"]
ALDIR = aln_opts["aligned_direc"]
NPROC = aln_opts["align_threads"]
SAMT_NPROC = aln_opts["samtools_threads"]
MIN_FRAG_LEN = aln_opts["min_fragment_length"]
MAX_FRAG_LEN = aln_opts["max_fragment_length"]
WRITE_UNAL = aln_opts["write_unaligned_reads_to_bam"]

SEQ_DB = conf_dict_global["genome"]["genome_base"]
# base directory for all code for ipod analysis
BINDIR = conf_dict_global["general"]["bindir"]
RAWDIR = conf_dict_global["general"]["rawdir"]
STARTDIR = os.getcwd()

# define some functions that will be used in the rest of the script
def run_bowtie(prefix, phredbase, db=SEQ_DB, pe=True):
    '''Run alignment using bowtie2

    Args:
    -----
    prefix : str
        Characters to prepend to files generated during alignment
    phredbase : int
        Quality score encoding
    db : str
        Path to bowtie2 database to serve as reference for alignment
    pe : bool
        Whether to run alignemtn on paired end reads

    Returns:
    --------
    None
        Runs alignment and generates sam file
    '''

    fwd = os.path.join(PROCDIR, prefix+F_READ_SUFFIX)
    rev = os.path.join(PROCDIR, prefix+R_READ_SUFFIX)
    fwd_unpaired = os.path.join( PROCDIR, prefix + F_UP_READ_SUFFIX )
    rev_unpaired = os.path.join( PROCDIR, prefix + R_UP_READ_SUFFIX )
    samout = os.path.join(ALDIR,prefix+"_bowtie2.sam")
    inspect_cmd = f"bowtie2-inspect {db} "\
        f"> bowtie2_inspect.log "\
        f"2> bowtie2_inspect.err"
    res = subprocess.call(inspect_cmd, shell=True)

    if res != 0:
        raise NoRefException(db)
  
    if pe:
        cmdline = f"bowtie2 -x {db} -1 {fwd} -2 {rev} "\
            f"-U {fwd_unpaired},{rev_unpaired} -S {samout} "\
            f"-q --end-to-end --very-sensitive "\
            f"-p {NPROC} --phred{phredbase} "\
            f"--fr -I {MIN_FRAG_LEN} -X {MAX_FRAG_LEN} --dovetail"
    else:
        cmdline = f"bowtie2 -x {db} "\
            f"-U {fwd} -S {samout} "\
            f"-q --end-to-end --very-sensitive "\
            f"-p {NPROC} --phred{phredbase}"

    if not WRITE_UNAL:
        cmdline += " --no-unal"

    # no need to call TMPDIR.cleanup() when used in a context manager like so
    with tempfile.TemporaryDirectory() as TMPDIR:
        cmdline += f" --temp-directory {TMPDIR} "\
            f"> {prefix}_bowtie2.log "\
            f"2> {prefix}_bowtie2.err"
        print("\n{}\n".format(cmdline))
        res = subprocess.call(cmdline, shell=True)

    if res == 0:
        print(f"samfile {samout} successfully generated")
    else:
        print(
            f"*** Encountered an error while running bowtie2. Check {prefix}_bowtie2.err"
        )

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
    cmdline1 = f"samtools view -bS {samname} > {bamname_un}"
    cmdline2 = f"samtools sort -@ {SAMT_NPROC} -o {bamname} {bamname_un}"
    cmdline3 = f"samtools index {bamname}"
    print("\n{}\n".format(cmdline1))
    print("\n{}\n".format(cmdline2))
    print("\n{}\n".format(cmdline3))
    retcode1 = subprocess.call(cmdline1,shell=True)
    retcode2 = subprocess.call(cmdline2,shell=True)
    retcode3 = subprocess.call(cmdline3,shell=True)
    if retcode1 == 0 and retcode2 == 0 and retcode3 == 0:

        print(f"Safe to remove samfile {prefix}")
        os.remove(samname)
        os.remove(bamname_un)
    else:
        sys.exit(f"*** Encountered an error while postprocessing {prefix}")

    # deduplicate alignments if umi handling is to be done
    if UMI:
        umi_sep = ":"
        if UMI_METHOD == "5-prime":
            umi_sep = "_"
        dedup_bamname = os.path.join(ALDIR, prefix+"_bowtie2_dedup.bam")
        dedup_pre = os.path.join(ALDIR, prefix+"_dedup_stats")
        dedup_cmd = f"umi_tools dedup --umi-separator {umi_sep} "\
            f"-I {bamname} --output-stats={dedup_pre} -S {dedup_bamname}"
        retcode4 = subprocess.call(dedup_cmd, shell=True)
        if retcode4 == 0:
            print(
                f"Deduplication ran without error, running:\n" \
                f"mv {dedup_bamname} {bamname}"
            )
            subprocess.call(f"mv {dedup_bamname} {bamname}", shell=True)
        else:
            sys.exit(f"*** Encountered an error while deduplicating {prefix}")

conf_dict = toml.load(sys.argv[1]) # this is the condition-level conf file
samp_types = conf_dict["general"]["sample_types"]

for samp_type in samp_types:

    print("Running alignments for {} sample.".format(samp_type))
    # gather this sample type's info from conf file
    samp_info = conf_dict[samp_type]
    freads = samp_info["R1_raw_files"]
    rreads = samp_info["R2_raw_files"]
    if "adapter_seqs" in samp_info:
        fwd_adapts = samp_info["adapter_seqs"]
        rev_adapts = samp_info["adapter_seqs"]
    if "R1_adapter_seqs" in samp_info:
        fwd_adapts = samp_info["R1_adapter_seqs"]
    if "R2_adapter_seqs" in samp_info:
        rev_adapts = samp_info["R2_adapter_seqs"]
    rep_names = samp_info["sample_names"]
    sample_direc = samp_info["directory"]

    # move into this sample's directory
    os.chdir(sample_direc)
    # set up the needed directories if they are not already present
    if not(os.path.isdir(ALDIR)):
        os.mkdir(ALDIR)

    # if the rawdir option is not "None" (note a string, not a None object)
    #   then determine whether it's already a symlink. If it's not a symlink
    #   already, create the symlink within the data directory
    if RAWDIR != "None":
        if not(os.path.islink('raw')):
            if not(os.path.isdir('raw')):
                os.symlink(conf_dict_global["general"]["rawdir"], "raw")

    # loop over each replicate's information and do preprocessing
    for i in range(len(rep_names)):
        samp_dict = {
            "ffile": freads[i],
            "rfile": rreads[i],
            "fwd_adapseq": fwd_adapts[i],
            "rev_adapseq": rev_adapts[i],
            "phredbase": conf_dict_global["general"]["phredbase"],
            "outprefix": rep_names[i],
            "pe": PE,
        }
        run_bowtie(
            samp_dict["outprefix"],
            samp_dict["phredbase"],
            db=SEQ_DB,
            pe=PE,
        )
        postprocess_bowtie(samp_dict["outprefix"])

    # when finished with this sample's alignments, move back up
    #  to condition-level directory
    os.chdir(STARTDIR)

